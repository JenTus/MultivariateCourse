Bivariate Correspondence Analysis
================

Bivariate Correspondence Analysis
---------------------------------

The data SMOKING1.txt contains a two-dimensional frequency table, where the employees of a company have been categorized according to their position (5 cat- egories: SM = Senior Managers, JM= Junior Managers, SE = Senior Employees, JE = Junior Employees, SC = Secretaries). The smoking of the employees have 4 categories (None, Light, Medium, Heavy).

Q1
--

Form the row and column profiles.

``` r
smok <- read.table("data/SMOKING.txt", header = TRUE, row.names = 1)
smok <- smok[-dim(smok)[1], -dim(smok)[2]] #remove the total col/row
# View(smok)
D <- as.matrix(smok)
prop.table(D, 1) # row profiles
```

    ##         None     Light    Medium      Heavy
    ## SM 0.3636364 0.1818182 0.2727273 0.18181818
    ## JM 0.2222222 0.1666667 0.3888889 0.22222222
    ## SE 0.4901961 0.1960784 0.2352941 0.07843137
    ## JE 0.2045455 0.2727273 0.3750000 0.14772727
    ## SC 0.4000000 0.2400000 0.2800000 0.08000000

``` r
prop.table(D, 2) # col profiles
```

    ##          None      Light    Medium Heavy
    ## SM 0.06557377 0.04444444 0.0483871  0.08
    ## JM 0.06557377 0.06666667 0.1129032  0.16
    ## SE 0.40983607 0.22222222 0.1935484  0.16
    ## JE 0.29508197 0.53333333 0.5322581  0.52
    ## SC 0.16393443 0.13333333 0.1129032  0.08

``` r
row.sum <- apply(smok, 1, sum)
row.profile <- sweep(smok, 1, row.sum, "/")

col.sum <- apply(smok, 2, sum)
col.profile <- sweep(smok, 2, col.sum, "/")
```

chi-square test
---------------

``` r
n <- sum(smok)
v1 <- matrix(colSums(smok), nrow = 1)
v2 <- matrix(rowSums(smok), ncol = 1)
E <- v2 %*% v1 / n # theoretical frequencies under independence
chi <- sum((D - E)^2 / E) # chi-square statistics
I <- dim(D)[1]
J <- dim(D)[2]
pchisq(chi, df = ((I-1)*(J-1)), lower.tail = F) # p-value
```

    ## [1] 0.1718348

``` r
chisq.test(smok)
```

    ## Warning in chisq.test(smok): Chi-squared approximation may be incorrect

    ## 
    ##  Pearson's Chi-squared test
    ## 
    ## data:  smok
    ## X-squared = 16.442, df = 12, p-value = 0.1718

Analysis ca
-----------

``` r
library(ca)
n <- sum(smok)
smok.ca <- ca(smok)
names(smok.ca)
```

    ##  [1] "sv"         "nd"         "rownames"   "rowmass"    "rowdist"   
    ##  [6] "rowinertia" "rowcoord"   "rowsup"     "colnames"   "colmass"   
    ## [11] "coldist"    "colinertia" "colcoord"   "colsup"     "N"         
    ## [16] "call"

``` r
smok.ca$sv # sqrt(lambda_i)
```

    ## [1] 0.27342111 0.10008587 0.02033652

``` r
smok.ca$rowdist # Row chi-square distances to centroid
```

    ## [1] 0.2165590 0.3569210 0.3807790 0.2400247 0.2161692

``` r
smok.ca$rownames
```

    ## [1] "SM" "JM" "SE" "JE" "SC"

``` r
smok.ca$coldist
```

    ## [1] 0.3944897 0.1739958 0.1981274 0.3551093

``` r
smok.ca$colnames
```

    ## [1] "None"   "Light"  "Medium" "Heavy"

``` r
smok.ca$rowcoord # scaled coordinates phi_ji/sqrt(lambda_i)
```

    ##          Dim1       Dim2       Dim3
    ## SM -0.2405388 -1.9357079  3.4903231
    ## JM  0.9471047 -2.4309584 -1.6573725
    ## SE -1.3919733 -0.1065076 -0.2535221
    ## JE  0.8519895  0.5769437  0.1625337
    ## SC -0.7354557  0.7884353 -0.3973677

``` r
smok.ca$colcoord
```

    ##              Dim1        Dim2        Dim3
    ## None   -1.4384714 -0.30465911 -0.04378737
    ## Light   0.3637463  1.40943267  1.08170100
    ## Medium  0.7180168  0.07352795 -1.26172451
    ## Heavy   1.0744451 -1.97595989  1.28885615

``` r
smok.ca$rowcoord[1,]*smok.ca$sv*1000 # coordinates in the summary
```

    ##       Dim1       Dim2       Dim3 
    ##  -65.76838 -193.73700   70.98103

``` r
sqrt(sum((smok.ca$rowcoord[1,]*smok.ca$sv)^2)) # distance from the center for SM
```

    ## [1] 0.216559

``` r
#innertia is the chi squared statistic divided by n
smok.ca$rowinertia 
```

    ## [1] 0.002672932 0.011881177 0.038314129 0.026268627 0.006052995

``` r
smok.ca$colinertia 
```

    ## [1] 0.049186258 0.007058828 0.012610242 0.016334533

``` r
sum(smok.ca$rowinertia)
```

    ## [1] 0.08518986

``` r
sum(smok.ca$colinertia)
```

    ## [1] 0.08518986

``` r
sum(smok.ca$sv^2)
```

    ## [1] 0.08518986

``` r
chi/n
```

    ## [1] 0.08518986

``` r
#take an example
sum((smok[1,] - E[1,])^2/E[1,])/n
```

    ## [1] 0.002672932

``` r
smok.ca$rowinertia[1]
```

    ## [1] 0.002672932

``` r
#proportional values, the same as in the summary
smok.ca$rowinertia/sum(smok.ca$sv^2) * 1000
```

    ## [1]  31.37618 139.46703 449.74987 308.35392  71.05300

``` r
# relative marginal frequencies of the original table
smok.ca$rowmass
```

    ## [1] 0.05699482 0.09326425 0.26424870 0.45595855 0.12953368

``` r
rowSums(smok/n)
```

    ##         SM         JM         SE         JE         SC 
    ## 0.05699482 0.09326425 0.26424870 0.45595855 0.12953368

``` r
margin.table(as.matrix(smok), 1)/sum(smok)
```

    ##         SM         JM         SE         JE         SC 
    ## 0.05699482 0.09326425 0.26424870 0.45595855 0.12953368

``` r
smok.ca$colmass
```

    ## [1] 0.3160622 0.2331606 0.3212435 0.1295337

``` r
colSums(smok/n)
```

    ##      None     Light    Medium     Heavy 
    ## 0.3160622 0.2331606 0.3212435 0.1295337

``` r
# the original table
smok.ca$N
```

    ##      [,1] [,2] [,3] [,4]
    ## [1,]    4    2    3    2
    ## [2,]    4    3    7    4
    ## [3,]   25   10   12    4
    ## [4,]   18   24   33   13
    ## [5,]   10    6    7    2

``` r
# ctr, the contribution in forming that ca-component, SM
margins <- smok.ca$rowmass
margins[1]*smok.ca$rowcoord[1,2]^2 # ctr
```

    ## [1] 0.2135576

``` r
# quality of representation --- Squared correlation
# Squared correlation of SM, quality of representation
d1 = (smok.ca$rowcoord[1,]*smok.ca$sv)^2 # sum(phi_ij^2)
d1/sum(d1)
```

    ##       Dim1       Dim2       Dim3 
    ## 0.09223203 0.80033639 0.10743158

Q2
--

How much of the variation is explained by the combination of components 1 and 3. Give the answer in percentages relative to the total variation.

``` r
library(ca)
smok.ca <- ca(smok)
smok.ca$sv^2/sum(smok.ca$sv^2)
```

    ## [1] 0.877558731 0.117586535 0.004854734

``` r
summary(smok.ca)
```

    ## 
    ## Principal inertias (eigenvalues):
    ## 
    ##  dim    value      %   cum%   scree plot               
    ##  1      0.074759  87.8  87.8  **********************   
    ##  2      0.010017  11.8  99.5  ***                      
    ##  3      0.000414   0.5 100.0                           
    ##         -------- -----                                 
    ##  Total: 0.085190 100.0                                 
    ## 
    ## 
    ## Rows:
    ##     name   mass  qlt  inr    k=1 cor ctr    k=2 cor ctr  
    ## 1 |   SM |   57  893   31 |  -66  92   3 | -194 800 214 |
    ## 2 |   JM |   93  991  139 |  259 526  84 | -243 465 551 |
    ## 3 |   SE |  264 1000  450 | -381 999 512 |  -11   1   3 |
    ## 4 |   JE |  456 1000  308 |  233 942 331 |   58  58 152 |
    ## 5 |   SC |  130  999   71 | -201 865  70 |   79 133  81 |
    ## 
    ## Columns:
    ##     name   mass  qlt  inr    k=1 cor ctr    k=2 cor ctr  
    ## 1 | None |  316 1000  577 | -393 994 654 |  -30   6  29 |
    ## 2 | Lght |  233  984   83 |   99 327  31 |  141 657 463 |
    ## 3 | Medm |  321  983  148 |  196 982 166 |    7   1   2 |
    ## 4 | Hevy |  130  995  192 |  294 684 150 | -198 310 506 |

Q3
--

Produce the BCA graph with respect to the first two components.

``` r
library
```

    ## function (package, help, pos = 2, lib.loc = NULL, character.only = FALSE, 
    ##     logical.return = FALSE, warn.conflicts = TRUE, quietly = FALSE, 
    ##     verbose = getOption("verbose")) 
    ## {
    ##     testRversion <- function(pkgInfo, pkgname, pkgpath) {
    ##         if (is.null(built <- pkgInfo$Built)) 
    ##             stop(gettextf("package %s has not been installed properly\n", 
    ##                 sQuote(pkgname)), call. = FALSE, domain = NA)
    ##         R_version_built_under <- as.numeric_version(built$R)
    ##         if (R_version_built_under < "3.0.0") 
    ##             stop(gettextf("package %s was built before R 3.0.0: please re-install it", 
    ##                 sQuote(pkgname)), call. = FALSE, domain = NA)
    ##         current <- getRversion()
    ##         if (length(Rdeps <- pkgInfo$Rdepends2)) {
    ##             for (dep in Rdeps) if (length(dep) > 1L) {
    ##                 target <- dep$version
    ##                 res <- if (is.character(target)) {
    ##                   do.call(dep$op, list(as.numeric(R.version[["svn rev"]]), 
    ##                     as.numeric(sub("^r", "", dep$version))))
    ##                 }
    ##                 else {
    ##                   do.call(dep$op, list(current, as.numeric_version(target)))
    ##                 }
    ##                 if (!res) 
    ##                   stop(gettextf("This is R %s, package %s needs %s %s", 
    ##                     current, sQuote(pkgname), dep$op, target), 
    ##                     call. = FALSE, domain = NA)
    ##             }
    ##         }
    ##         if (R_version_built_under > current) 
    ##             warning(gettextf("package %s was built under R version %s", 
    ##                 sQuote(pkgname), as.character(built$R)), call. = FALSE, 
    ##                 domain = NA)
    ##         platform <- built$Platform
    ##         r_arch <- .Platform$r_arch
    ##         if (.Platform$OS.type == "unix") {
    ##             if (!nzchar(r_arch) && length(grep("\\w", platform)) && 
    ##                 !testPlatformEquivalence(platform, R.version$platform)) 
    ##                 stop(gettextf("package %s was built for %s", 
    ##                   sQuote(pkgname), platform), call. = FALSE, 
    ##                   domain = NA)
    ##         }
    ##         else {
    ##             if (nzchar(platform) && !grepl("mingw", platform)) 
    ##                 stop(gettextf("package %s was built for %s", 
    ##                   sQuote(pkgname), platform), call. = FALSE, 
    ##                   domain = NA)
    ##         }
    ##         if (nzchar(r_arch) && file.exists(file.path(pkgpath, 
    ##             "libs")) && !file.exists(file.path(pkgpath, "libs", 
    ##             r_arch))) 
    ##             stop(gettextf("package %s is not installed for 'arch = %s'", 
    ##                 sQuote(pkgname), r_arch), call. = FALSE, domain = NA)
    ##     }
    ##     checkLicense <- function(pkg, pkgInfo, pkgPath) {
    ##         L <- tools:::analyze_license(pkgInfo$DESCRIPTION["License"])
    ##         if (!L$is_empty && !L$is_verified) {
    ##             site_file <- path.expand(file.path(R.home("etc"), 
    ##                 "licensed.site"))
    ##             if (file.exists(site_file) && pkg %in% readLines(site_file)) 
    ##                 return()
    ##             personal_file <- path.expand("~/.R/licensed")
    ##             if (file.exists(personal_file)) {
    ##                 agreed <- readLines(personal_file)
    ##                 if (pkg %in% agreed) 
    ##                   return()
    ##             }
    ##             else agreed <- character()
    ##             if (!interactive()) 
    ##                 stop(gettextf("package %s has a license that you need to accept in an interactive session", 
    ##                   sQuote(pkg)), domain = NA)
    ##             lfiles <- file.path(pkgpath, c("LICENSE", "LICENCE"))
    ##             lfiles <- lfiles[file.exists(lfiles)]
    ##             if (length(lfiles)) {
    ##                 message(gettextf("package %s has a license that you need to accept after viewing", 
    ##                   sQuote(pkg)), domain = NA)
    ##                 readline("press RETURN to view license")
    ##                 encoding <- pkgInfo$DESCRIPTION["Encoding"]
    ##                 if (is.na(encoding)) 
    ##                   encoding <- ""
    ##                 if (encoding == "latin1") 
    ##                   encoding <- "cp1252"
    ##                 file.show(lfiles[1L], encoding = encoding)
    ##             }
    ##             else {
    ##                 message(gettextf("package %s has a license that you need to accept:\naccording to the DESCRIPTION file it is", 
    ##                   sQuote(pkg)), domain = NA)
    ##                 message(pkgInfo$DESCRIPTION["License"], domain = NA)
    ##             }
    ##             choice <- menu(c("accept", "decline"), title = paste("License for", 
    ##                 sQuote(pkg)))
    ##             if (choice != 1) 
    ##                 stop(gettextf("license for package %s not accepted", 
    ##                   sQuote(package)), domain = NA, call. = FALSE)
    ##             dir.create(dirname(personal_file), showWarnings = FALSE)
    ##             writeLines(c(agreed, pkg), personal_file)
    ##         }
    ##     }
    ##     checkNoGenerics <- function(env, pkg) {
    ##         nenv <- env
    ##         ns <- .getNamespace(as.name(pkg))
    ##         if (!is.null(ns)) 
    ##             nenv <- asNamespace(ns)
    ##         if (exists(".noGenerics", envir = nenv, inherits = FALSE)) 
    ##             TRUE
    ##         else {
    ##             length(objects(env, pattern = "^\\.__T", all.names = TRUE)) == 
    ##                 0L
    ##         }
    ##     }
    ##     checkConflicts <- function(package, pkgname, pkgpath, nogenerics, 
    ##         env) {
    ##         dont.mind <- c("last.dump", "last.warning", ".Last.value", 
    ##             ".Random.seed", ".Last.lib", ".onDetach", ".packageName", 
    ##             ".noGenerics", ".required", ".no_S3_generics", ".Depends", 
    ##             ".requireCachedGenerics")
    ##         sp <- search()
    ##         lib.pos <- match(pkgname, sp)
    ##         ob <- objects(lib.pos, all.names = TRUE)
    ##         if (!nogenerics) {
    ##             these <- ob[substr(ob, 1L, 6L) == ".__T__"]
    ##             gen <- gsub(".__T__(.*):([^:]+)", "\\1", these)
    ##             from <- gsub(".__T__(.*):([^:]+)", "\\2", these)
    ##             gen <- gen[from != package]
    ##             ob <- ob[!(ob %in% gen)]
    ##         }
    ##         fst <- TRUE
    ##         ipos <- seq_along(sp)[-c(lib.pos, match(c("Autoloads", 
    ##             "CheckExEnv"), sp, 0L))]
    ##         for (i in ipos) {
    ##             obj.same <- match(objects(i, all.names = TRUE), ob, 
    ##                 nomatch = 0L)
    ##             if (any(obj.same > 0)) {
    ##                 same <- ob[obj.same]
    ##                 same <- same[!(same %in% dont.mind)]
    ##                 Classobjs <- grep("^\\.__", same)
    ##                 if (length(Classobjs)) 
    ##                   same <- same[-Classobjs]
    ##                 same.isFn <- function(where) vapply(same, exists, 
    ##                   NA, where = where, mode = "function", inherits = FALSE)
    ##                 same <- same[same.isFn(i) == same.isFn(lib.pos)]
    ##                 not.Ident <- function(ch, TRAFO = identity, ...) vapply(ch, 
    ##                   function(.) !identical(TRAFO(get(., i)), TRAFO(get(., 
    ##                     lib.pos)), ...), NA)
    ##                 if (length(same)) 
    ##                   same <- same[not.Ident(same)]
    ##                 if (length(same) && identical(sp[i], "package:base")) 
    ##                   same <- same[not.Ident(same, ignore.environment = TRUE)]
    ##                 if (length(same)) {
    ##                   if (fst) {
    ##                     fst <- FALSE
    ##                     packageStartupMessage(gettextf("\nAttaching package: %s\n", 
    ##                       sQuote(package)), domain = NA)
    ##                   }
    ##                   msg <- .maskedMsg(same, pkg = sQuote(sp[i]), 
    ##                     by = i < lib.pos)
    ##                   packageStartupMessage(msg, domain = NA)
    ##                 }
    ##             }
    ##         }
    ##     }
    ##     if (verbose && quietly) 
    ##         message("'verbose' and 'quietly' are both true; being verbose then ..")
    ##     if (!missing(package)) {
    ##         if (is.null(lib.loc)) 
    ##             lib.loc <- .libPaths()
    ##         lib.loc <- lib.loc[dir.exists(lib.loc)]
    ##         if (!character.only) 
    ##             package <- as.character(substitute(package))
    ##         if (length(package) != 1L) 
    ##             stop("'package' must be of length 1")
    ##         if (is.na(package) || (package == "")) 
    ##             stop("invalid package name")
    ##         pkgname <- paste("package", package, sep = ":")
    ##         newpackage <- is.na(match(pkgname, search()))
    ##         if (newpackage) {
    ##             pkgpath <- find.package(package, lib.loc, quiet = TRUE, 
    ##                 verbose = verbose)
    ##             if (length(pkgpath) == 0L) {
    ##                 txt <- if (length(lib.loc)) 
    ##                   gettextf("there is no package called %s", sQuote(package))
    ##                 else gettext("no library trees found in 'lib.loc'")
    ##                 if (logical.return) {
    ##                   warning(txt, domain = NA)
    ##                   return(FALSE)
    ##                 }
    ##                 else stop(txt, domain = NA)
    ##             }
    ##             which.lib.loc <- normalizePath(dirname(pkgpath), 
    ##                 "/", TRUE)
    ##             pfile <- system.file("Meta", "package.rds", package = package, 
    ##                 lib.loc = which.lib.loc)
    ##             if (!nzchar(pfile)) 
    ##                 stop(gettextf("%s is not a valid installed package", 
    ##                   sQuote(package)), domain = NA)
    ##             pkgInfo <- readRDS(pfile)
    ##             testRversion(pkgInfo, package, pkgpath)
    ##             if (!package %in% c("datasets", "grDevices", "graphics", 
    ##                 "methods", "splines", "stats", "stats4", "tcltk", 
    ##                 "tools", "utils") && isTRUE(getOption("checkPackageLicense", 
    ##                 FALSE))) 
    ##                 checkLicense(package, pkgInfo, pkgpath)
    ##             if (is.character(pos)) {
    ##                 npos <- match(pos, search())
    ##                 if (is.na(npos)) {
    ##                   warning(gettextf("%s not found on search path, using pos = 2", 
    ##                     sQuote(pos)), domain = NA)
    ##                   pos <- 2
    ##                 }
    ##                 else pos <- npos
    ##             }
    ##             .getRequiredPackages2(pkgInfo, quietly = quietly)
    ##             deps <- unique(names(pkgInfo$Depends))
    ##             if (packageHasNamespace(package, which.lib.loc)) {
    ##                 if (isNamespaceLoaded(package)) {
    ##                   newversion <- as.numeric_version(pkgInfo$DESCRIPTION["Version"])
    ##                   oldversion <- as.numeric_version(getNamespaceVersion(package))
    ##                   if (newversion != oldversion) {
    ##                     res <- try(unloadNamespace(package))
    ##                     if (inherits(res, "try-error")) 
    ##                       stop(gettextf("Package %s version %s cannot be unloaded", 
    ##                         sQuote(package), oldversion, domain = "R-base"))
    ##                   }
    ##                 }
    ##                 tt <- try({
    ##                   ns <- loadNamespace(package, c(which.lib.loc, 
    ##                     lib.loc))
    ##                   env <- attachNamespace(ns, pos = pos, deps)
    ##                 })
    ##                 if (inherits(tt, "try-error")) 
    ##                   if (logical.return) 
    ##                     return(FALSE)
    ##                   else stop(gettextf("package or namespace load failed for %s", 
    ##                     sQuote(package)), call. = FALSE, domain = NA)
    ##                 else {
    ##                   on.exit(detach(pos = pos))
    ##                   nogenerics <- !.isMethodsDispatchOn() || checkNoGenerics(env, 
    ##                     package)
    ##                   if (warn.conflicts && !exists(".conflicts.OK", 
    ##                     envir = env, inherits = FALSE)) 
    ##                     checkConflicts(package, pkgname, pkgpath, 
    ##                       nogenerics, ns)
    ##                   on.exit()
    ##                   if (logical.return) 
    ##                     return(TRUE)
    ##                   else return(invisible(.packages()))
    ##                 }
    ##             }
    ##             else stop(gettextf("package %s does not have a namespace and should be re-installed", 
    ##                 sQuote(package)), domain = NA)
    ##         }
    ##         if (verbose && !newpackage) 
    ##             warning(gettextf("package %s already present in search()", 
    ##                 sQuote(package)), domain = NA)
    ##     }
    ##     else if (!missing(help)) {
    ##         if (!character.only) 
    ##             help <- as.character(substitute(help))
    ##         pkgName <- help[1L]
    ##         pkgPath <- find.package(pkgName, lib.loc, verbose = verbose)
    ##         docFiles <- c(file.path(pkgPath, "Meta", "package.rds"), 
    ##             file.path(pkgPath, "INDEX"))
    ##         if (file.exists(vignetteIndexRDS <- file.path(pkgPath, 
    ##             "Meta", "vignette.rds"))) 
    ##             docFiles <- c(docFiles, vignetteIndexRDS)
    ##         pkgInfo <- vector("list", 3L)
    ##         readDocFile <- function(f) {
    ##             if (basename(f) %in% "package.rds") {
    ##                 txt <- readRDS(f)$DESCRIPTION
    ##                 if ("Encoding" %in% names(txt)) {
    ##                   to <- if (Sys.getlocale("LC_CTYPE") == "C") 
    ##                     "ASCII//TRANSLIT"
    ##                   else ""
    ##                   tmp <- try(iconv(txt, from = txt["Encoding"], 
    ##                     to = to))
    ##                   if (!inherits(tmp, "try-error")) 
    ##                     txt <- tmp
    ##                   else warning("'DESCRIPTION' has an 'Encoding' field and re-encoding is not possible", 
    ##                     call. = FALSE)
    ##                 }
    ##                 nm <- paste0(names(txt), ":")
    ##                 formatDL(nm, txt, indent = max(nchar(nm, "w")) + 
    ##                   3)
    ##             }
    ##             else if (basename(f) %in% "vignette.rds") {
    ##                 txt <- readRDS(f)
    ##                 if (is.data.frame(txt) && nrow(txt)) 
    ##                   cbind(basename(gsub("\\.[[:alpha:]]+$", "", 
    ##                     txt$File)), paste(txt$Title, paste0(rep.int("(source", 
    ##                     NROW(txt)), ifelse(nzchar(txt$PDF), ", pdf", 
    ##                     ""), ")")))
    ##                 else NULL
    ##             }
    ##             else readLines(f)
    ##         }
    ##         for (i in which(file.exists(docFiles))) pkgInfo[[i]] <- readDocFile(docFiles[i])
    ##         y <- list(name = pkgName, path = pkgPath, info = pkgInfo)
    ##         class(y) <- "packageInfo"
    ##         return(y)
    ##     }
    ##     else {
    ##         if (is.null(lib.loc)) 
    ##             lib.loc <- .libPaths()
    ##         db <- matrix(character(), nrow = 0L, ncol = 3L)
    ##         nopkgs <- character()
    ##         for (lib in lib.loc) {
    ##             a <- .packages(all.available = TRUE, lib.loc = lib)
    ##             for (i in sort(a)) {
    ##                 file <- system.file("Meta", "package.rds", package = i, 
    ##                   lib.loc = lib)
    ##                 title <- if (nzchar(file)) {
    ##                   txt <- readRDS(file)
    ##                   if (is.list(txt)) 
    ##                     txt <- txt$DESCRIPTION
    ##                   if ("Encoding" %in% names(txt)) {
    ##                     to <- if (Sys.getlocale("LC_CTYPE") == "C") 
    ##                       "ASCII//TRANSLIT"
    ##                     else ""
    ##                     tmp <- try(iconv(txt, txt["Encoding"], to, 
    ##                       "?"))
    ##                     if (!inherits(tmp, "try-error")) 
    ##                       txt <- tmp
    ##                     else warning("'DESCRIPTION' has an 'Encoding' field and re-encoding is not possible", 
    ##                       call. = FALSE)
    ##                   }
    ##                   txt["Title"]
    ##                 }
    ##                 else NA
    ##                 if (is.na(title)) 
    ##                   title <- " ** No title available ** "
    ##                 db <- rbind(db, cbind(i, lib, title))
    ##             }
    ##             if (length(a) == 0L) 
    ##                 nopkgs <- c(nopkgs, lib)
    ##         }
    ##         dimnames(db) <- list(NULL, c("Package", "LibPath", "Title"))
    ##         if (length(nopkgs) && !missing(lib.loc)) {
    ##             pkglist <- paste(sQuote(nopkgs), collapse = ", ")
    ##             msg <- sprintf(ngettext(length(nopkgs), "library %s contains no packages", 
    ##                 "libraries %s contain no packages"), pkglist)
    ##             warning(msg, domain = NA)
    ##         }
    ##         y <- list(header = NULL, results = db, footer = NULL)
    ##         class(y) <- "libraryIQR"
    ##         return(y)
    ##     }
    ##     if (logical.return) 
    ##         TRUE
    ##     else invisible(.packages())
    ## }
    ## <bytecode: 0x1e3ec60>
    ## <environment: namespace:base>

``` r
plot(smok.ca, arrows=c(TRUE, T), map = "symmetric")
```

![](home6_files/figure-markdown_github/unnamed-chunk-5-1.png)

``` r
plot(smok.ca, arrows=c(TRUE, T), map = "symmetric", dim = c(1,3))
```

![](home6_files/figure-markdown_github/unnamed-chunk-5-2.png)

``` r
plot(smok.ca, arrows = c(T, T), map = "rowprincipal")
```

![](home6_files/figure-markdown_github/unnamed-chunk-5-3.png)

Q4
--

How much of the variation of the modality Heavy is explained by the first two components? How much of the variation of the modality Medium is explained by the first two components? Give the answers in percentages relative to the total variation. Hint: use the quality of representation.

``` r
# the variation of the modality Heavy explained by the first two components
h <- (smok.ca$colcoord[4,]*smok.ca$sv)^2
sum(h[c(1,2)])/sum(h)
```

    ## [1] 0.994552

``` r
# the variation of the modality Medium explained by the first two components
m <- (smok.ca$colcoord[3,]*smok.ca$sv)^2
sum(m[c(1,2)])/sum(m)
```

    ## [1] 0.9832277
