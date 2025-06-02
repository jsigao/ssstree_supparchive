## Modified BEAST Program
We provide an executable for a modified version of the `BEAST` program, which implements some additional functionalities, including the recording of site likelihood values at each sampled iteration and some extensions of the `TreeSummary` tool.
The modified source code (from which we compiled the executable) is located in [the `master` branch, commit `6c4d255`](https://github.com/jsigao/beast-mcmc/commit/6c4d255829c78a293e2ce3c1729c56939d64a09a) of a forked `BEAST` source-code repository.
Briefly, this executable can be run by invoking:
```
java -jar beast.jar -working ./analyses/beast/ebov_dud17/run01.xml
```
See [`BEAST` official tutorial website](http://beast.community/index.html) for further details about running `BEAST` analyses using a `Java` executable.

Note that the XML scripts we provide in the `analyses` subdirectory can be ran directly with the released version of `BEAST` v.1.10.5 as the site-likelihood logging sections are commented out in those scripts to ensure compatibility.
The `get_sitell_posthoc.R` script we provide in the `scripts` subdirectory can then be used to compute the site-likelihood values in a post-hoc manner.