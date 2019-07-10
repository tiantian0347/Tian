#include "cache.h"
#include "repository.h"
#include "refs.h"
#include "object.h"
#include "commit.h"
#include "tag.h"
#include "packfile.h"
#include "object-store.h"

/*
 * Create the file "path" by writing to a temporary file and renaming
 * it into place. The contents of the file come from "generate", which
 * should return non-zero if it encounters an error.
 */
static int update_info_file(char *path, int (*generate)(FILE *))
{
	char *tmp = mkpathdup("%s_XXXXXX", path);
	int ret = -1;
	int fd = -1;
	FILE *fp = NULL, *to_close;

	safe_create_leading_directories(path);
	fd = git_mkstemp_mode(tmp, 0666);
	if (fd < 0)
		goto out;
	to_close = fp = fdopen(fd, "w");
	if (!fp)
		goto out;
	fd = -1;
	ret = generate(fp);
	if (ret)
		goto out;
	fp = NULL;
	if (fclose(to_close))
		goto out;
	if (adjust_shared_perm(tmp) < 0)
		goto out;
	if (rename(tmp, path) < 0)
		goto out;
	ret = 0;

out:
	if (ret) {
		error_errno("unable to update %s", path);
		if (fp)
			fclose(fp);
		else if (fd >= 0)
			close(fd);
		unlink(tmp);
	}
	free(tmp);
	return ret;
}

static int add_info_ref(const char *path, const struct object_id *oid,
			int flag, void *cb_data)
{
	FILE *fp = cb_data;
	struct object *o = parse_object(the_repository, oid);
	if (!o)
		return -1;

	if (fprintf(fp, "%s	%s\n", oid_to_hex(oid), path) < 0)
		return -1;

	if (o->type == OBJ_TAG) {
		o = deref_tag(the_repository, o, path, 0);
		if (o)
			if (fprintf(fp, "%s	%s^{}\n",
				oid_to_hex(&o->oid), path) < 0)
				return -1;
	}
	return 0;
}

static int generate_info_refs(FILE *fp)
{
	return for_each_ref(add_info_ref, fp);
}

static int update_info_refs(void)
{
	char *path = git_pathdup("info/refs");
	int ret = update_info_file(path, generate_info_refs);
	free(path);
	return ret;
}

/* packs */
static struct pack_info {
	struct packed_git *p;
	int old_num;
	int new_num;
} **info;
static int num_pack;

static struct pack_info *find_pack_by_name(const char *name)
{
	int i;
	for (i = 0; i < num_pack; i++) {
		struct packed_git *p = info[i]->p;
		if (!strcmp(pack_basename(p), name))
			return info[i];
	}
	return NULL;
}

/* Returns non-zero when we detect that the info in the
 * old file is useless.
 */
static int parse_pack_def(const char *packname, int old_cnt)
{
	struct pack_info *i = find_pack_by_name(packname);
	if (i) {
		i->old_num = old_cnt;
		return 0;
	}
	else {
		/* The file describes a pack that is no longer here */
		return 1;
	}
}

/* Returns non-zero when we detect that the info in the
 * old file is useless.
 */
static int read_pack_info_file(const char *infofile)
{
	FILE *fp;
	struct strbuf line = STRBUF_INIT;
	int old_cnt = 0;
	int stale = 1;

	fp = fopen_or_warn(infofile, "r");
	if (!fp)
		return 1; /* nonexistent is not an error. */

	while (strbuf_getline(&line, fp) != EOF) {
		const char *arg;

		if (!line.len)
			continue;

		if (skip_prefix(line.buf, "P ", &arg)) {
			/* P name */
			if (parse_pack_def(arg, old_cnt++))
				goto out_stale;
		} else if (line.buf[0] == 'D') {
			/* we used to emit D but that was misguided. */
			goto out_stale;
		} else if (line.buf[0] == 'T') {
			/* we used to emit T but nobody uses it. */
			goto out_stale;
		} else {
			error("unrecognized: %s", line.buf);
		}
	}
	stale = 0;

 out_stale:
	strbuf_release(&line);
	fclose(fp);
	return stale;
}

static int compare_info(const void *a_, const void *b_)
{
	struct pack_info *const *a = a_;
	struct pack_info *const *b = b_;

	if (0 <= (*a)->old_num && 0 <= (*b)->old_num)
		/* Keep the order in the original */
		return (*a)->old_num - (*b)->old_num;
	else if (0 <= (*a)->old_num)
		/* Only A existed in the original so B is obviously newer */
		return -1;
	else if (0 <= (*b)->old_num)
		/* The other way around. */
		return 1;

	/* then it does not matter but at least keep the comparison stable */
	if ((*a)->p == (*b)->p)
		return 0;
	else if ((*a)->p < (*b)->p)
		return -1;
	else
		return 1;
}

static void init_pack_info(const char *infofile, int force)
{
	struct packed_git *p;
	int stale;
	int i = 0;

	for (p = get_all_packs(the_repository); p; p = p->next) {
		/* we ignore things on alternate path since they are
		 * not available to the pullers in general.
		 */
		if (!p->pack_local)
			continue;
		i++;
	}
	num_pack = i;
	info = xcalloc(num_pack, sizeof(struct pack_info *));
	for (i = 0, p = get_all_packs(the_repository); p; p = p->next) {
		if (!p->pack_local)
			continue;
		assert(i < num_pack);
		info[i] = xcalloc(1, sizeof(struct pack_info));
		info[i]->p = p;
		info[i]->old_num = -1;
		i++;
	}

	if (infofile && !force)
		stale = read_pack_info_file(infofile);
	else
		stale = 1;

	for (i = 0; i < num_pack; i++)
		if (stale)
			info[i]->old_num = -1;

	/* renumber them */
	QSORT(info, num_pack, compare_info);
	for (i = 0; i < num_pack; i++)
		info[i]->new_num = i;
}

static void free_pack_info(void)
{
	int i;
	for (i = 0; i < num_pack; i++)
		free(info[i]);
	free(info);
}

static int write_pack_info_file(FILE *fp)
{
	int i;
	for (i = 0; i < num_pack; i++) {
		if (fprintf(fp, "P %s\n", pack_basename(info[i]->p)) < 0)
			return -1;
	}
	if (fputc('\n', fp) == EOF)
		return -1;
	return 0;
}

static int update_info_packs(int force)
{
	char *infofile = mkpathdup("%s/info/packs", get_object_directory());
	int ret;

	init_pack_info(infofile, force);
	ret = update_info_file(infofile, write_pack_info_file);
	free_pack_info();
	free(infofile);
	return ret;
}

/* public */
int update_server_info(int force)
{
	/* We would add more dumb-server support files later,
	 * including index of available pack files and their
	 * intended audiences.
	 */
	int errs = 0;

	errs = errs | update_info_refs();
	errs = errs | update_info_packs(force);

	/* remove leftover rev-cache file if there is any */
	unlink_or_warn(git_path("info/rev-cache"));

	return errs;
}
