	
	/*** test for bernoulli proportions, best-arm selection with δ-PAC guarantees ***/

	var bernoulli_pac = function(delta_value) {

		var stoppingTime;
		var delta = delta_value; // the error guarantee we want
		var x_data = [];
		var y_data = [];
		var n = 0;
		var S_x = 0;
		var S_y = 0;
		var finished = false;
		var L_an;

		/** public functions **/

		this.getResults = function() {
			var L_an = LikH0(S_x, S_y, n);
			return {
				'S_x' : S_x,
				'S_y' : S_y,
				'L_an' : L_an,
				'finished' : finished,
				'n' : n
			};
		}

		// get confidence interval (only when test is done)
		this.confInterval = function(samples) {
			if (!finished) {
				return undefined;
			}
			if (!samples) samples = 10000;

			// get unbiased result
			var ests = this.estimate();

			var outcomes = [];
			// simulate n outcomes
			for (var i = 0;i < samples;i++) {
				var res = simulateResult(ests[0],ests[1]);
				var time = res[3];
				outcomes[i] = [res[1]/time, res[2]/time];
			}
			outcomes.sort(function(a,b){return (a[0]-a[1])-(b[0]-b[1]);})

			// bias corrected bootstrap confidence interval
			var outcomes_diff = [];
			var lower_count = 0;
			for (var i = 0;i < outcomes.length;i++) {
				outcomes_diff[i] = outcomes[i][0] - outcomes[i][1];
				if (outcomes_diff[i] < ((S_x/n)-(S_y/n))) lower_count += 1;
			}
			//console.log("lower count:"+lower_count)
			var b = jStat.jStat.normal.inv(lower_count/samples,0,1);
			//console.log(b);
			var upper_n = Math.floor((samples+1)*jStat.jStat.normal.cdf(2*b + 1.96,0,1));
			var lower_n = Math.floor((samples+1)*jStat.jStat.normal.cdf(2*b - 1.96,0,1));
			//console.log("lower_n:"+lower_n)
			//console.log("upper_n:"+upper_n)
			var lower_est = outcomes[lower_n];
			var upper_est = outcomes[upper_n];

			// bias correct the lower and upper estimates
			var lower_est_bc = optimize2d(lower_est, biasFun(), lower_est, 0.005, 16400, 590000, 0.02, 0, 1, false);
			var upper_est_bc = optimize2d(upper_est, biasFun(), upper_est, 0.005, 16400, 590000, 0.02, 0, 1, false);

			return [(lower_est_bc[0]-lower_est_bc[1]),(upper_est_bc[0]-upper_est_bc[1])];
		}

		// get estimate (only when test is done)
		  // use bias-reduction
		this.estimate = function() {
			if (!finished) {
				return undefined;
			}
			var ests = optimize2d([S_x/n, S_y/n], biasFun(), [S_x/n, S_y/n], 0.005, 16400, 590000, 0.02, 0, 1, false);
			// TODO : should we include std.dev.?
			return [ests[0], ests[1], ests[0]-ests[1]];
		}
		
		// get sequence of data
		this.getData = function() {
			return [x_data, y_data];
		}
		
		// add single or paired datapoint (control or treatment)
			// returns true if test is finished
		this.addData = function(points) {
			if (finished) {
				if (typeof points[0] === 'number') x_data.push(points[0]);
				if (typeof points[1] === 'number') y_data.push(points[1]);
			} else {
				if (typeof points[0] === 'number' && typeof points[1] === 'number') {
					if (x_data.length == y_data.length) {
						S_x += points[0];
						S_y += points[1];
						n += 1;
					} else if (x_data.length > y_data.length) {
						S_y += points[1];
						S_x += x_data[n];
						n += 1;
					} else {
						S_x += points[0];
						S_y += y_data[n];
						n += 1;
					}
					x_data.push(points[0])
					y_data.push(points[1])
				} else if (typeof points[0] === 'number') {
					if (x_data.length < y_data.length) {
						S_x += points[0];
						S_y += y_data[n];
						n += 1;
					}
					x_data.push(points[0]);
				} else if (typeof points[1] === 'number') {
					if (x_data.length > y_data.length) {
						S_y += points[1];
						S_x += x_data[n];
						n += 1;
					}
					y_data.push(points[1]);
				}
			}
			
			var result = checkTest(S_x, S_y, n);
			if (result) {
				finished = true;
				stoppingTime = n;
				return result;
			}
		}

		// get expected samplesize for some parameters
		/*this.expectedSamplesize = function(p1, p2, samples) {
			// simulate it enough times
			if (!samples) samples = 10000;
			console.log("calculating expected samplesize via simulation");
			var times = [];
			for (var i = 0;i < samples;i++) {
				var res = simulateResult(p1,p2)
				times.push(res[3]);
			}
			return mean(times);
		}*/

		this.expectedSamplesize = function(p1, p2, samples) {
			// simulate it enough times
			if (!samples) samples = 10000;
			console.log("calculating expected samplesize via simulation");
			var times = [];
			for (var i = 0;i < samples;i++) {
				var res = simulateResult(p1,p2)
				times.push(res[3]);
			}
			times.sort(function(a,b){return a-b});
			return [times[samples*0.05], mean(times), times[samples*0.95]];
		}

		/** private functions **/
		
		var biasFun = function() {
			var outfun = function(pt, n) {
				var results_p1 = []
				var results_p2 = []
				for (var i = 0;i < n;i++) {
					// generate sequences
					var res = simulateResult(pt[0], pt[1]);
					results_p1.push( res[1]/res[3] );
					results_p2.push( res[2]/res[3] );
				}
				return [results_p1, results_p2];
			}
			return outfun;
		}
		
		var checkTest = function(S_x, S_y, n) {
			// check if test should be stopped
			var L_an = LikH0(S_x, S_y, n);
			if (L_an >= Math.log(2*n)/delta) {
				if (S_x > S_y) {
					return 'X';
				} else {
					return 'Y';
				}
			}
			return undefined
		}

		var LikH0 = functions['bernoulli_pac']['l_an'];

		var generate = function(p) {
			if (Math.random() < p) {return 1;} else {return 0;}
		}
		
		var simulateResult = function(p1, p2) {
			var finished = false;
			var time = 0;
			var S_x = 0;
			var S_y = 0;
			var result;
			while (!finished) {
				S_x += generate(p1);
				S_y += generate(p2);
				time += 1;
				// test it
				var result = checkTest(S_x, S_y, time);
				if (result) finished = true;
			}
			return [result, S_x, S_y, time];
		}

		// get test variables
		this.properties = {
			'delta' : delta,
		}
	}

	// private functions

	var bernoulli_pac_LR_H0 = function(S_x, S_y, n) {
		var equal_mle = (S_x+S_y)/(2*n);
		// calculate unconstrained MLE, i.e. p1 and p2 can be unequal 
		var unc_mle_x = S_x/n;
		var unc_mle_y = S_y/n;

		var likRatio = Math.exp( (logOp(S_x,unc_mle_x) + logOp(n-S_x,1-unc_mle_x) + logOp(S_y,unc_mle_y) + logOp(n-S_y,1-unc_mle_y)) - (logOp(S_x,equal_mle) + logOp(n-S_x,1-equal_mle) + logOp(S_y,equal_mle) + logOp(n-S_y,1-equal_mle)));
		return likRatio;
	}

	functions['bernoulli_pac'] = {
		'l_an' : bernoulli_pac_LR_H0, // this is in bernoulli.js
	}
	
	tests['bernoulli_pac'] = bernoulli_pac;