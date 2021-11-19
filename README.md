# SeqAn Lib Template [![build status][1]][2] [![codecov][3]][4]
<!--
    Above uses reference-style links with numbers.
    See also https://github.com/adam-p/markdown-here/wiki/Markdown-Cheatsheet#links.

    For example, `[![build status][1]][2]` evaluates to the following:
        `[link_text][2]`
        `[2]` is a reference to a link, i.e. `[link_text](https://...)`

        `[link_text]` = `[![build status][1]]`
        `[1]` is once again a reference to a link - this time an image, i.e. `[![build status](https://...)]
        `![build status]` is the text that should be displayed if the linked resource (`[1]`) is not available

    `[![build status][1]][2]` hence means:
    Show the picture linked under `[1]`. In case it cannot be displayed, show the text "build status" instead.
    The picture, or alternative text, should link to `[2]`.
-->

<!--
    This is the CI badge image:
        `https://img.shields.io/github/workflow/status/` - we do not use GitHub's badges as they are not customisable.
        `/rrahn/lib_template/` - owner/repository
        `CI%20on%20Linux` - name of the workflow as encoded URL (e.g., whitespace = %20)
        `master` - branch to show
        `?style=flat&logo=github` - use a GitHub-style badge
        `&label=lib_template%20CI` - text on the badge
        `"Open GitHub actions page"` - this text will be shown on hover
-->
[1]: https://img.shields.io/github/workflow/status/rrahn/lib_template/CI%20on%20Linux/master?style=flat&logo=github&label=lib_template%20CI "Open GitHub actions page"
<!--
    This is the CI badge link:
        `https://github.com/rrahn/lib_template/actions` - actions page of owner(seqan)/repository(lib_template)
        `?query=branch%3Amaster` - only show actions that ran on the mater branch
-->
[2]: https://github.com/rrahn/lib_template/actions?query=branch%3Amaster
<!--
    This is the Codecov badge image:
        Codecov offers badges: https://app.codecov.io/gh/rrahn/lib_template/settings/badge
        While being logged in into Codecov, navigate to Settings->Badge and copy the markdown badge.
        Copy the image part of the markdown badge here.
    `"Open Codecov page"` - this text will be shown on hover
-->
[3]: https://codecov.io/gh/rrahn/lib_template/branch/master/graph/badge.svg?token=V82JRCXF0K "Open Codecov page"
<!--
    This is the Codecov badge link:
        Codecov offers badges: https://app.codecov.io/gh/rrahn/lib_template/settings/badge
        While being logged in into Codecov, navigate to Settings->Badge and copy the markdown badge.
        Copy the URL part of the markdown badge here.
-->
[4]: https://codecov.io/gh/rrahn/lib_template

This is a library template for SeqAn developers contributing new algorithms or data structures to the SeqAn project.
You can easily clone this repository and modify the existing code to your needs.
It provides the elementary set-up for all SeqAn3 libraries in order to provide a consistent infrastructure.

After cloning this repository, please rename all instances of `lib_template` or `LIB_TEMPLATE` with the respective project name `<project_name>` in all files and also update corresponding file names.
Note that for the SeqAn project we only accept header-only files.

Instructions:
1. clone this repository: `git clone --recurse-submodules https://github.com/rrahn/lib_template.git lib_template`
2. rename all `lib_template` occurrences in all files and also all filenames with the `<project_name>`.
3. create a build directory and visit it: `mkdir build && cd build`
4. configure unit tests cmake: `cmake ../<project_name>/test/unit`
5. build the unit tests: `make`
6. run the unit tests: `ctest`

Library structure:
1. build_system: Defines cmake export target and manages installation of the library.
2. doc: Contains the infrastructure necessary to build the API documentation.
3. include: Contains the header files for your project.
4. submodule: External git submodules this library depends on.
5. test: The test directory containing a separate module for unit tests, benchmark tests, coverage tests, test data, header tests, and a cmake directory containing additional cmake files.
