
############## host on github ################
library(GetoptLong)

### create repos

username = "simplifyEnrichment"
token = 


setwd("/icgc/dkfzlsdf/analysis/B080/guz/simplifyGO_test/examples")

library(GetoptLong)
files = list.files()
dir = files[file.info(files)$isdir]

for(repo in dir) {
	# -XDELETE for delete repos
	cmd = qq("curl -u '@{username}:@{token}' https://api.github.com/user/repos -d \"{\\\"name\\\":\\\"@{repo}\\\"}\"")
	qqcat("creating repo for @{repo}\n")
	system(cmd)
}

##### push
for(repo in dir) {

	qqcat("------------ upload for @{repo} -------------\n")

	setwd(qq("/icgc/dkfzlsdf/analysis/B080/guz/simplifyGO_test/examples/@{repo}"))
	if(!file.exists(".git")) {
		system("git init")
		cat("rds\n", file = ".gitignore")
		qqcat("figures for @{repo}", file = "readme.md")
		system("git add readme.md")
		system("git commit -m 'add readme.md'")
		system("git branch gh-pages")
	}

	system("git add --all")
	system(qq("git commit -m 'add files for @{repo}'"))
	system(qq("git push https://@{username}:@{token}@github.com/simplifyEnrichment/@{repo}.git master"))
	system("git checkout gh-pages")
	system("git merge master")
	system(qq("git push https://@{username}:@{token}@github.com/simplifyEnrichment/@{repo}.git gh-pages"))
	system("git checkout master")
}

################### similarity comparisons ###################
setwd("/icgc/dkfzlsdf/analysis/B080/guz/simplifyGO_test/compare_similarity")


library(GetoptLong)
files = list.files()
dir = files[file.info(files)$isdir]
dir = setdiff(dir, "rds")

for(repo in dir) {
	# -XDELETE for delete repos
	cmd = qq("curl -u '@{username}:@{token}' https://api.github.com/user/repos -d \"{\\\"name\\\":\\\"cmp_sim_@{repo}\\\"}\"")
	qqcat("creating repo for cmp_sim_@{repo}\n")
	system(cmd)
}

##### push
for(repo in dir) {

	qqcat("------------ upload for cmp_sim_@{repo} -------------\n")

	setwd(qq("/icgc/dkfzlsdf/analysis/B080/guz/simplifyGO_test/compare_similarity/@{repo}"))
	if(!file.exists(".git")) {
		system("git init")
		cat("rds\n", file = ".gitignore")
		qqcat("figures for cmp_sim_@{repo}", file = "readme.md")
		system("git add readme.md")
		system("git commit -m 'add readme.md'")
		system("git branch gh-pages")
	}

	system("git add --all")
	system(qq("git commit -m 'add files for cmp_sim_@{repo}'"))
	system(qq("git push https://@{username}:@{token}@github.com/simplifyEnrichment/cmp_sim_@{repo}.git master"))
	system("git checkout gh-pages")
	system("git merge master")
	system(qq("git push https://@{username}:@{token}@github.com/simplifyEnrichment/cmp_sim_@{repo}.git gh-pages"))
	system("git checkout master")
}
