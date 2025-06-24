<!-- omit in toc -->
# Contributing to seabreeze

First off, thanks for taking the time to contribute! 

All types of contributions are encouraged and valued. See the [Table of Contents](#table-of-contents) for different ways to help and details about how this project handles them. Please make sure to read the relevant section before making your contribution. It will make it a lot easier for us maintainers and smooth out the experience for all involved. The community looks forward to your contributions. 

> And if you like the project, but just don't have time to contribute, that's fine. There are other easy ways to support the project and show your appreciation, which we would also be very happy about:
> - Star the project
> - Tweet about it
> - Refer this project in your project's readme
> - Mention the project at local meetups and tell your friends/colleagues

<!-- omit in toc -->
## Table of Contents

- [I Have a Question](#i-have-a-question)
- [I Want To Contribute](#i-want-to-contribute)
- [Reporting Bugs](#reporting-bugs)
- [Suggesting Enhancements](#suggesting-enhancements)
- [Improving The Documentation](#improving-the-documentation)



## I Have a Question

> If you want to ask a question, we assume that you have read the available [Documentation](https://barricklab.github.io/seabreeze/).

Before you ask a question, it is best to search for existing [Issues](https://github.com/barricklab/seabreeze/issues) that might help you. In case you have found a suitable issue and still need clarification, you can write your question in this issue. It is also advisable to search the internet for answers first.

If you then still feel the need to ask a question and need clarification, we recommend the following:

- Open an [Issue](https://github.com/barricklab/seabreeze/issues/new).
- Provide as much context as you can about what you're running into.
- Provide project and platform versions. In particular, the version of `conda` being used. 

We will then take care of the issue as soon as possible.


## I Want To Contribute

> ### Legal Notice <!-- omit in toc -->
> When contributing to this project, you must agree that you have authored 100% of the content, that you have the necessary rights to the content and that the content you contribute may be provided under the project licence.

We love pull requests from everyone. To get started, first fork, then clone the repo. After introducing any changes, please run the workflow test. We use `pytest` for testing, have implemented two types of tests:

- unit tests: These are for individual scripts. Build the test environment first

```
conda env create --file test/unit_tests/environment.yml  --name seabreeze_test
```

Then activate the environment, and run the unit tests. Note that the unit tests should be run from the seabreeze root directory. 

```
conda activate seabreeze_test
pytest test/unit_tests/*.py
```

`ISEScan`, one of the dependencies, is not compatible with some ARM-based systems. 

- workflow tests: These test the CLI itself, and ensures that it correctly calls and executes snakemake workflows. These tests can be run from the seabreeze execution environment itself.

```
conda activate seabreeze
pytest test/*.py
```



 To test that your contribution is compatible with the snakemake pipeline, use the example data in the `example/` folder. Push to your fork and submit a pull request.

### Reporting Bugs

<!-- omit in toc -->
#### Before Submitting a Bug Report

A good bug report shouldn't leave others needing to chase you up for more information. Therefore, we ask you to investigate carefully, collect information and describe the issue in detail in your report. Please complete the following steps in advance to help us fix any potential bug as fast as possible.

- Make sure that you are using the latest version.
- Determine if your bug is really a bug and not an error on your side e.g. using incompatible environment components/versions (Make sure that you have read the [documentation](https://barricklab.github.io/seabreeze/). If you are looking for support, you might want to check [this section](#i-have-a-question)).
- To see if other users have experienced (and potentially already solved) the same issue you are having, check if there is not already a bug report existing for your bug or error in the [bug tracker](https://github.com/barricklab/seabreeze/issues?q=label%3Abug).

- Collect information about the bug:
  - Stack trace (Traceback). Snakemake generates log files in the hidden `.snakemake` directory. Each rule in the workflow also generates a log file, in `data/logs`. If the workflow is failing at a particular rule, we recommend including the log file for that rule as well. 
  - OS, Platform and Version (Windows, Linux, macOS, x86, ARM)
  - Possibly your input and the output
  - Can you reliably reproduce the issue? 

<!-- omit in toc -->
#### How Do I Submit a Good Bug Report?

> You must never report security related issues, vulnerabilities or bugs including sensitive information to the issue tracker, or elsewhere in public. Instead sensitive bugs must be sent by email to <ira.zibbu@gmail.com>.
<!-- You may add a PGP key to allow the messages to be sent encrypted as well. -->

We use GitHub issues to track bugs and errors. If you run into an issue with the project:

- Open an [Issue](https://github.com/barricklab/seabreeze/issues/new). (Since we can't be sure at this point whether it is a bug or not, we ask you not to talk about a bug yet and not to label the issue.)
- Explain the behavior you would expect and the actual behavior.
- Please provide as much context as possible and describe the *reproduction steps* that someone else can follow to recreate the issue on their own. This usually includes your code. For good bug reports you should isolate the problem and create a reduced test case.
- Provide the information you collected in the previous section.

Once it's filed, we will work on addressing the issue.

<!-- You might want to create an issue template for bugs and errors that can be used as a guide and that defines the structure of the information to be included. If you do so, reference it here in the description. -->


### Suggesting Enhancements

This section guides you through submitting an enhancement suggestion for seabreeze, **including completely new features and minor improvements to existing functionality**. Following these guidelines will help maintainers and the community to understand your suggestion and find related suggestions.

<!-- omit in toc -->
#### Before Submitting an Enhancement

- Make sure that you are using the latest version.
- Read the [documentation](https://barricklab.github.io/seabreeze/) carefully and find out if the functionality is already covered, maybe by an individual configuration.
- Perform a [search](https://github.com/barricklab/seabreeze/issues) to see if the enhancement has already been suggested. If it has, add a comment to the existing issue instead of opening a new one.

<!-- omit in toc -->
#### How Do I Submit a Good Enhancement Suggestion?

Enhancement suggestions are tracked as [GitHub issues](https://github.com/barricklab/seabreeze/issues).

- Use a **clear and descriptive title** for the issue to identify the suggestion.
- Provide a **step-by-step description of the suggested enhancement** in as many details as possible.
- **Describe the current behavior** and **explain which behavior you expected to see instead** and why. At this point you can also tell which alternatives do not work for you.
- You may want to **include screenshots or screen recordings** which help you demonstrate the steps or point out the part which the suggestion is related to. You can use [LICEcap](https://www.cockos.com/licecap/) to record GIFs on macOS and Windows, and the built-in [screen recorder in GNOME](https://help.gnome.org/users/gnome-help/stable/screen-shot-record.html.en) or [SimpleScreenRecorder](https://github.com/MaartenBaert/ssr) on Linux. <!-- this should only be included if the project has a GUI -->
- **Explain why this enhancement would be useful** to most seabreeze users. You may also want to point out the other projects that solved it better and which could serve as inspiration.

## Improving the Documentation

If you believe the documentation is incomplete or could use additional information to guide other users, we highly welcome your contribution. seabreeze uses [Mkdocs](https://www.mkdocs.org/) to generate documentation from Markdown files that exist in the `documentation\` directory. To contribute, fork the repository, make the desired changes in the `documentation` directory and submit a pull request. 

<!-- omit in toc -->
## Attribution
This guide is based on the [contributing.md](https://contributing.md/generator)!
