SeGLiR
======

**SeGLiR** is a javascript library for rapid A/B-testing with **Sequential Generalized Likelihood Ratio Tests**.

Sequential generalized likelihood ratio tests is a family of [sequential hypothesis tests](http://en.wikipedia.org/wiki/Sequential_analysis), i.e. tests that stop as soon as a significant result has been detected in accordance with some stopping rule. Compared to classical fixed samplesize tests, sequential generalized likelihood ratio tests may give a significant decrease in the needed samplesize, while keeping the same type-1 and type-2 error guarantees. Below is a comparison of the expected samplesize for fixed samplesize and sequential glr tests at the same levels. For a more detailed explanation of this family of tests, see [the reference](http://auduno.github.io/seglir/documentation.html#math).

![Expected samplesize comparison](https://dl.dropboxusercontent.com/u/10557805/samplesize2b.png)

SeGLiR currently contain these tests:
* comparing two bernoulli proportions (for instance for comparing conversion on a website)
* comparing two normal means (with equal, known or unknown variance)
* choosing the best arm in a multi-armed bandit setting (with δ-PAC guarantees)

Improvements and additions are welcome, take a look at [issues](https://github.com/auduno/seglir/issues) for outstanding issues.

### Installation ###

The recommended way is to install SeGLiR via node package manager (install node.js if you don't already have it):

```
npm install seglir
```

### Usage ###

If you're not running SeGLiR in node, note that SeGLiR requires [*jStat.js*](https://github.com/jstat/jstat) available.

```javascript
var glr = require(seglir);
// create an instance of a two-sided test comparing bernoulli proportions, with indifference region with size 0.01, alpha-level = 0.05, beta-level = 0.10
var test = new glr.test("bernoulli", "two-sided", 0.01, 0.05, 0.10);
// add data as it comes in
...
test.addData({x : 0});
test.addData({x : 0, y : 0});
test.addData({y : 1});
...
// when function returns string 'true' or 'false', the test is concluded
test.getResults()
// get bias-adjusted estimates
test.estimates()
```

For a complete function reference, see the [reference](http://auduno.github.io/seglir/documentation.html).

### Building from source ###

Make sure you have [grunt](http://gruntjs.com/) and [node](http://nodejs.org/download/) installed.
To install the development dependencies run ```npm install``` and to build it run ```grunt``` in the root directory.

### License ###

SeGLiR is distributed under the [MIT License](http://www.opensource.org/licenses/MIT).
