
############## host on github ################
library(GetoptLong)
library(gh)

username = "simplifyEnrichment"
token = 

delete_repo = function(repo, dir, where = c("remote", "local")) {
	qqcat("delete repo @{repo}\n")
	if("remote" %in% where) {
		gh("DELETE /repos/:owner/:repo", owner = username, repo = repo, .token = token)
	}

	if("local" %in% where) {
		owd = getwd()
		on.exit(setwd(owd))
		setwd(dir)
		if(file.exists(".git")) {
			unlink(".git", recursive = TRUE, force = TRUE)
		}
	}
}

create_repo = function(repo, dir, ignore = c("rds", "rds/", "^.*", "Rplots.pdf"),
	where = c("remote", "local")) {
	
	all_repos = gh(qq("GET /users/@{username}/repos"))
	all_repos = vapply(all_repos, "[[", "", "name")
	if(repo %in% all_repos) {
		delete_repo(repo, dir)
	} else {
		delete_repo(repo, dir, where = "local")
	}

	qqcat("creating repo for @{repo}\n")
	if("remote" %in% where) {
		new_repo = gh("POST /user/repos", name = repo, owner = username, .token = token)
	}

	if("local" %in% where) {
		owd = getwd()
		on.exit(setwd(owd))
		setwd(dir)

		system("git init")
		cat(ignore, sep = "\n", file = ".gitignore")
		system("git add .gitignore")
		system("git commit -m 'add .gitignore'")

		if(!file.exists("readme.md")) {
			qqcat("figures for @{repo}", file = "readme.md")
			system("git add readme.md")
			system("git commit -m 'add readme.md'")
		}

		branch = system("git branch", intern = TRUE)
		if(!any(grepl("gh-pages", branch))) {
			system("git branch gh-pages")
		}
	}
}

update_repo = function(repo, dir) {
	qqcat("------------ upload for @{repo} -------------\n")

	owd = getwd()
	on.exit(setwd(owd))
	setwd(dir)

	system("git add --all")
	system(qq("git commit -m 'add files for @{repo}'"))
	system(qq("git push https://@{username}:@{token}@github.com/@{username}/@{repo}.git master"))
	system("git checkout gh-pages")
	system("git merge master")
	system(qq("git push https://@{username}:@{token}@github.com/@{username}/@{repo}.git gh-pages"))
	system("git checkout master")
}

sub_dirs = function(dir, pattern = NULL) {
	owd = getwd()
	on.exit(setwd(owd))
	setwd(dir)

	files = list.files(pattern = pattern)
	files[file.info(files)$isdir]
}

dir = sub_dirs("/icgc/dkfzlsdf/analysis/B080/guz/simplifyGO_test/examples")

for(repo in dir) {
	create_repo(repo, qq("/icgc/dkfzlsdf/analysis/B080/guz/simplifyGO_test/examples/@{repo}"))
}

##### push
for(repo in dir) {
	update_repo(repo, qq("/icgc/dkfzlsdf/analysis/B080/guz/simplifyGO_test/examples/@{repo}"))
}

################### similarity comparisons ###################
dir = sub_dirs("/icgc/dkfzlsdf/analysis/B080/guz/simplifyGO_test/compare_similarity")
dir = setdiff(dir, "rds")

for(repo in dir) {
	create_repo(qq("cmp_sim_@{repo}"), qq("/icgc/dkfzlsdf/analysis/B080/guz/simplifyGO_test/compare_similarity/@{repo}"))
}

##### push
for(repo in dir) {
	update_repo(qq("cmp_sim_@{repo}"), qq("/icgc/dkfzlsdf/analysis/B080/guz/simplifyGO_test/compare_similarity/@{repo}"))
}

################# test partition methods ###########
create_repo("test_partition_methods")
update_repo("test_partition_methods", "/icgc/dkfzlsdf/analysis/B080/guz/simplifyGO_test/test_partition_methods")


create_repo("compare_similarity", "/icgc/dkfzlsdf/analysis/B080/guz/simplifyGO_test/compare_similarity",
	ignore = sub_dirs("/icgc/dkfzlsdf/analysis/B080/guz/simplifyGO_test/compare_similarity", pattern = "EBI|random|rds"))
update_repo("compare_similarity", "/icgc/dkfzlsdf/analysis/B080/guz/simplifyGO_test/compare_similarity")


create_repo("examples", "/icgc/dkfzlsdf/analysis/B080/guz/simplifyGO_test/examples",
	ignore = sub_dirs("/icgc/dkfzlsdf/analysis/B080/guz/simplifyGO_test/examples"))
update_repo("examples", "/icgc/dkfzlsdf/analysis/B080/guz/simplifyGO_test/examples")
