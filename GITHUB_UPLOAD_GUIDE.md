# Upload and use BayesianDesignAna

The repository is named `ST1_1_Bayesian_Design_Ana`, while the R package itself
is named `BayesianDesignAna`. They intentionally have different names.

## Upload through the GitHub website (no Terminal required)

1. Sign in to https://github.com/JinhaoZou.
2. If you already have the `BayesianDesignAna` repository and want to rename it:
   open it, select **Settings**, change **Repository name** to
   `ST1_1_Bayesian_Design_Ana`, and choose **Rename**. Otherwise, create a
   repository with that name at https://github.com/new.
3. On your Mac, extract the package ZIP. Enter its **BayesianDesignAna** folder.
   Your main view should show `DESCRIPTION`, `NAMESPACE`, `R/`, `man/`, etc.
4. Review the `DESCRIPTION` maintainer email; replace `maintainer@example.org`
   with your email before publishing. Check ownership, licensing, and source
   content before making the repository public.
5. In Finder press **Command + Shift + .** to show dotfiles, including `.github`.
   Select and drag **all contents** of the extracted package folder onto the
   GitHub repository page under **Add file → Upload files**. Do **not** upload
   only the ZIP, and do not create an extra enclosing package folder.
6. Confirm that `DESCRIPTION` and `R/` appear in the repository's root.
   Enter a message like `Release v0.6.0` and commit the upload.
7. If you are updating an existing repository, review differences and replace
   changed files. Remove superseded files as needed; simply uploading new
   files will not remove old ones. Consider making the repository private until
   the package and statistical methods have been tested.

## Local installation on macOS

Extract the ZIP. In RStudio Console, specify the path to the folder containing
`DESCRIPTION` (adjust the path if yours differs):

```r
install.packages('remotes')
remotes::install_local('~/Downloads/BayesianDesignAna', force=TRUE,
                      upgrade='never')
library(BayesianDesignAna)
packageVersion('BayesianDesignAna')
```

If an older version is already loaded, use **Session → Restart R** before
installing and again afterward. You do not need to call `remove.packages()`
for each update.

## Installing from GitHub

After the repository is public (or with appropriate private-repository access):

```r
install.packages('remotes')
remotes::install_github('JinhaoZou/ST1_1_Bayesian_Design_Ana',
                        upgrade='never')
library(BayesianDesignAna)
```

## Quality checks

From the unzipped package's folder:

```r
install.packages('devtools')
devtools::test('~/Downloads/BayesianDesignAna')
devtools::check('~/Downloads/BayesianDesignAna')
```

Package checks do not demonstrate statistical validity. Compare model results
against reference analyses before scientific deployment.
