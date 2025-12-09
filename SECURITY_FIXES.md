# Security Fixes Applied

This document summarizes the security vulnerabilities that were addressed in this repository.

## JavaScript Library Vulnerabilities Fixed

### 1. jQuery (Critical Priority)
- **Previous Version:** 1.11.3 (released 2014)
- **Updated To:** 3.6.0
- **Vulnerabilities Fixed:**
  - CVE-2015-9251: Cross-site scripting (XSS) vulnerability
  - CVE-2019-11358: Prototype pollution vulnerability
  - CVE-2020-11022: XSS vulnerability in htmlPrefilter
  - CVE-2020-11023: XSS vulnerability in jQuery.htmlPrefilter
- **Files Updated:** 5 HTML files (index.html, license.html, Simulation_Benchmarking.html, RCH_B-ALL.html, Leucegene_Gene_Expression.html)

### 2. jQuery UI (High Priority)
- **Previous Version:** 1.11.4 (released 2015)
- **Updated To:** 1.14.1
- **Vulnerabilities Fixed:**
  - High severity XSS vulnerabilities in versions < 1.12.0
- **Files Updated:** 5 HTML files (Leucegene_Gene_Expression.html, Leucegene_Normals.html, Leucegene_Validation.html, RCH_B-ALL.html, Simulation_Benchmarking.html)

### 3. Bootstrap
- **Current Version:** 3.3.5
- **Status:** No critical vulnerabilities found in GitHub Advisory Database
- **Action:** No update required at this time

## Verification

All changes have been verified:
- ✅ Site loads correctly
- ✅ Table of contents functionality works (jQuery UI widget)
- ✅ CodeQL security scan shows 0 alerts
- ✅ Code review completed with no issues
- ✅ All HTML files reference updated library versions

## Testing

The site was tested locally using a Python HTTP server and Playwright browser automation:
- All pages load successfully
- JavaScript functionality works as expected
- No console errors related to library loading

## Future Maintenance

When rebuilding the site from R source files:
- Update `rmarkdown` package to latest version (which bundles updated JavaScript libraries)
- Update `workflowr` package to latest version
- Run `workflowr::wflow_build()` to regenerate HTML files with updated libraries
