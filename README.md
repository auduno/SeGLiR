SeGLiR
======

**SeGLiR** is a javascript library for rapid A/B-testing with **Sequential Generalized Likelihood Ratio Tests**.

For a brief explanation of this family of tests, see [this post](http://auduno.com/seglir/). SeGLiR currently contain tests for comparing two bernoulli proportions (for instance for comparing  conversion), but it's also planned to contain tests for comparing two normal means, as well as tests for choosing the best arm in a multi-armed bandit setting with δ-PAC guarantees (ref. Emilie Kaufman).

There are a number of improvements that can be done to this library, especially relating to the optimization algorithm for estimating a function minimum in two dimensions, so if you feel like helping out, take a look at the (issues)[].

### Usage ###

Download via node etc. Note that SeGLiR requires 'jStat'.

```javascript
var glr = require(seglir);
var test = new glr.test("bernoulli","two-sided",0.01,0.05,0.10);
// collect data
...
test.add([1,undefined]);
test.add([0,0]);
test.add([undefined,1]);
...
// when function returns 'true', 'false' or 'undefined' test is concluded
test.getResults()
// get bias-adjusted estimates 
test.estimates()
```

For a complete reference, see -..-

### Building from source ###

Make sure you have [grunt](http://gruntjs.com/) and [node](http://nodejs.org/download/) installed.
To install the development dependencies run ```npm install``` and to build it run ```grunt``` in the root directory.

### License ###

SeGLiR is distributed under the [MIT License](http://www.opensource.org/licenses/MIT).
