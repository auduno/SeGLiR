	
	/*** test for bernoulli proportions, best-arm selection with δ-PAC guarantees ***/

	var bernoulli_pac = function(delta_value) {

		// check input
		if (typeof(delta_value) != 'number' || delta_value <= 0 || delta_value >= 1) {
			console.log("parameter 'delta_value' must be a number between 0 and 1, input was : "+delta_value);
			return;
		}

		var delta = delta_value; // the error guarantee we want
		var x_data = [];
		var y_data = [];
		var n_x = 0;
		var n_y = 0;
		var S_x = 0;
		var S_y = 0;
		var finished = false;
		var L_an;

		/** public functions **/

		this.getResults = function() {
			var L_an = LikH0(S_x, S_y, n_x, n_y);
			return {
				'S_x' : S_x,
				'S_y' : S_y,
				'L_an' : L_an,
				'finished' : finished,
				'n_x' : n_x,
				'n_y' : n_y
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
				if (outcomes_diff[i] < ((S_x/n_x)-(S_y/n_y))) lower_count += 1;
			}
			//console.log("lower count:"+lower_count)
			var b = jStat.normal.inv(lower_count/samples,0,1);
			//console.log(b);
			var upper_n = Math.floor((samples+1)*jStat.normal.cdf(2*b + 1.96,0,1));
			var lower_n = Math.floor((samples+1)*jStat.normal.cdf(2*b - 1.96,0,1));
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
		this.estimate = function(max_samples) {
			if (!finished) {
				return undefined;
			}
			if (typeof(max_samples) == "undefined") {
				max_samples = 1500000;
			}
			var ests = optimize2d([S_x/n_x, S_y/n_y], biasFun(), [S_x/n_x, S_y/n_y], 0.005, 16400, max_samples, 0.02, 0, 1, true);
			// TODO : should we include std.dev.?
			return [ests[0], ests[1], ests[0]-ests[1]];
		}
		
		// get sequence of data
		this.getData = function() {
			return [x_data, y_data];
		}
		
		// add single or paired datapoint (control or treatment)
		this.addData = function(points) {
			var test = false;
			if (finished) {
				if (typeof points['x'] === 'number') x_data.push(points['x']);
				if (typeof points['y'] === 'number') y_data.push(points['y']);
			} else {
				if (typeof points['x'] === 'number' && typeof points['y'] === 'number') {
					if (x_data.length == y_data.length) {
						S_x += points['x'];
						S_y += points['y'];
					} else if (x_data.length > y_data.length) {
						S_y += points['y'];
						if (x_data.length == y_data.length+1) {
							S_x += points['x'];
						} else {
							S_x += x_data[n_x];
						}
					} else {
						S_x += points['x'];
						if (x_data.length+1 == y_data.length) {
							S_y += points['y'];
						} else {
							S_y += y_data[n_y];
						}
					}
					n_x += 1;
					n_y += 1;
					test = true;
					x_data.push(points['x'])
					y_data.push(points['y'])
				} else if (typeof points['x'] === 'number') {
					if (x_data.length == y_data.length) {
						S_x += points['x'];
						test = true;
						n_x += 1;
					} else if (x_data.length < y_data.length) {
						S_x += points['x'];
						test = true;
						n_x += 1;
						if (x_data.length+1 != y_data.length) {
							S_y += y_data[n_y];
							n_y += 1;
						}
					}
					x_data.push(points['x']);
				} else if (typeof points['y'] === 'number') {
					if (x_data.length == y_data.length) {
						S_y += points['y'];
						test = true;
						n_y += 1;
					} else if (x_data.length > y_data.length) {
						S_y += points['y'];
						test = true;
						n_y += 1;
						if (x_data.length != y_data.length+1) {
							S_x += x_data[n_x];
							n_x += 1;
						}
					} 
					y_data.push(points['y']);
				}
			}
			
			if (test) {
				var result = checkTest(S_x, S_y, n_x, n_y);
				if (result) {
					finished = true;
					return result;
				}
			}
		}

		// get expected samplesize for some parameters
		this.expectedSamplesize = function(p1, p2, samples) {
			// simulate it enough times
			if (!samples) samples = 10000;
			console.log("calculating expected samplesize via simulation");
			var times = [];
			for (var i = 0;i < samples;i++) {
				var res = simulateResult(p1,p2)
				times.push(res[3]);
			}
			return mean(times);
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
		
		var checkTest = function(S_x, S_y, n_x, n_y) {
			// check if test should be stopped
			var L_an = LikH0(S_x, S_y, n_x, n_y);
			if (L_an >= (Math.log(n_x + n_y)+1)/delta) {
				if (S_x/n_x > S_y/n_y) {
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
				var result = checkTest(S_x, S_y, time, time);
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

	var bernoulli_pac_LR_H0 = function(S_x, S_y, n_x, n_y) {
		var equal_mle = (S_x+S_y)/(n_x + n_y);
		// calculate unconstrained MLE, i.e. p1 and p2 can be unequal 
		var unc_mle_x = S_x/n_x;
		var unc_mle_y = S_y/n_y;

		var likRatio = Math.exp( (logOp(S_x,unc_mle_x) + logOp(n_x-S_x,1-unc_mle_x) + logOp(S_y,unc_mle_y) + logOp(n_y-S_y,1-unc_mle_y)) - (logOp(S_x,equal_mle) + logOp(n_x-S_x,1-equal_mle) + logOp(S_y,equal_mle) + logOp(n_y-S_y,1-equal_mle)));
		return likRatio;
	}

	functions['bernoulli_pac'] = {
		'l_an' : bernoulli_pac_LR_H0, // this is in bernoulli.js
	}
	
	tests['bernoulli_pac'] = bernoulli_pac;