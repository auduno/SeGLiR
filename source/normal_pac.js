	
	/*** test for comparing normal means with equal known variance, best-arm selection with δ-PAC guarantees ***/

	var normal_pac = function(delta_value, variance) {

		// check input
		if (typeof(delta_value) != 'number' || delta_value <= 0 || delta_value >= 1) {
			console.log("parameter 'delta_value' must be a number between 0 and 1, input was : "+delta_value);
			return;
		}
		if (typeof(variance) != 'number' || variance <= 0) {
			console.log("parameter 'variance' must be a number larger than 0, input was : "+variance);
			return;
		}

		var delta = delta_value;
		var var_value = variance;
		var x_data = [];
		var y_data = [];
		var n_x = 0;
		var n_y = 0;
		var S_x = 0;
		var S_y = 0;
		var S_x2 = 0;
		var S_y2 = 0;
		var finished = false;
		var L_an;

		/** public functions **/

		this.getResults = function() {
			var L_an = LikH0(S_x, S_y, S_x2, S_y2, n_x, n_y);
			return {
				'S_x' : S_x,
				'S_y' : S_y,
				'S_x2' : S_x2,
				'S_y2' : S_y2,
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
				var time = res[5];
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
			var ests = optimize2d([S_x/n_x, S_y/n_y], biasFun(), [S_x/n_x, S_y/n_y], 0.005, 16400, 590000, 0.02, 0, 1, false);
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
				if (typeof points[0] === 'number') x_data.push(points[0]);
				if (typeof points[1] === 'number') y_data.push(points[1]);
			} else {
				if (typeof points[0] === 'number' && typeof points[1] === 'number') {
					if (x_data.length == y_data.length) {
						S_x += points[0];
						S_y += points[1];
					} else if (x_data.length > y_data.length) {
						S_y += points[1];
						if (x_data.length == y_data.length+1) {
							S_x += points[0];
						} else {
							S_x += x_data[n_x];
						}
					} else {
						S_x += points[0];
						if (x_data.length+1 == y_data.length) {
							S_y += points[1];
						} else {
							S_y += y_data[n_y];
						}
					}
					n_x += 1;
					n_y += 1;
					test = true;
					x_data.push(points[0])
					y_data.push(points[1])
				} else if (typeof points[0] === 'number') {
					if (x_data.length == y_data.length) {
						S_x += points[0];
						test = true;
						n_x += 1;
					} else if (x_data.length < y_data.length) {
						S_x += points[0];
						test = true;
						n_x += 1;
						if (x_data.length+1 != y_data.length) {
							S_y += y_data[n_y];
							n_y += 1;
						}
					}
					x_data.push(points[0]);
				} else if (typeof points[1] === 'number') {
					if (x_data.length == y_data.length) {
						S_y += points[1];
						test = true;
						n_y += 1;
					} else if (x_data.length > y_data.length) {
						S_y += points[1];
						test = true;
						n_y += 1;
						if (x_data.length != y_data.length+1) {
							S_x += x_data[n_x];
							n_x += 1;
						}
					} 
					y_data.push(points[1]);
				}
			}
			
			if (test) {
				var result = checkTest(S_x, S_y, S_x2, S_y2, n_x, n_y);
				if (result) {
					finished = true;
					return result;
				}
			}
		}

		// get expected samplesize for some parameters
		this.expectedSamplesize = function(params_1, params_2, samples) {
			// simulate it enough times
			if (!samples) samples = 10000;
			console.log("calculating expected samplesize via simulation");
			var times = [];
			for (var i = 0;i < samples;i++) {
				var res = simulateResult(params_1,params_2)
				times.push(res[5]);
			}
			return mean(times);
		}
				
		/** private functions **/

		var generate = function(params) {
			var mean = params[0];
			var std = Math.sqrt(params[1]);
			return jStat.jStat.normal.sample(mean,std);
		}

		var simulateResult = function(params_1, params_2) {
			var sample1, sample2, result;
			var finished = false;
			var time = 0;
			var S_x = 0;
			var S_y = 0;
			var S_x2 = 0;
			var S_y2 = 0;
			while (!finished) {
				sample1 = generate(params_1);
				sample2 = generate(params_2);
				S_x += sample1;
				S_y += sample2;
				S_x2 += sample1*sample1;
				S_y2 += sample2*sample2;
				time += 1;
				// test it
				var result = checkTest(S_x, S_y, S_x2, S_y2, time, time);
				if (result) finished = true;
			}
			// return result, S_x, S_y, stoppingTime
			return [result[0], S_x, S_y, S_x2, S_y2, time, result[1]];
		}
		this.simulateResult = simulateResult;

		var checkTest = function(S_x, S_y, S_x2, S_y2, n_x, n_y) {
			// check if test should be stopped
			
			var L_an = LikH0(S_x, S_y, S_x2, S_y2, n_x, n_y, var_value);
			// better threshold
			var threshold = Math.pow((n_x+n_y+1)/(2*delta),(n_x+n_y+1)/(n_x+n_y));
			// more conservative threshold ? 
			//var threshold = (1/delta)*Math.pow(-Math.log(delta),3/4)*Math.pow(1+Math.log((n_x+n_y)/2),3/2)
			if (L_an >= threshold) {
				if (S_x/n_x > S_y/n_y) {
					return 'X';
				} else {
					return 'Y';
				}
			}
			return undefined
		}
		
		var biasFun = function() {
			var outfun = function(parameters, n) {
				var results_p1 = []
				var results_p2 = []
				for (var i = 0;i < n;i++) {
					// generate sequences
					var res = simulateResult(parameters[0], parameters[1]);
					results_p1.push( res[1]/res[5] );
					results_p2.push( res[2]/res[5] );
				}
				return [results_p1, results_p2];
			}
			return outfun;
		}

		var LikH0 = functions['normal_pac']['l_an']; // TODO : correct ?

		// get test variables
		this.properties = {
			'delta' : delta,
			'variance' : var_value
		}

		/*this.checkErrors = function(mu_1, mu_2, samples) {
			var errs = 0;
			if (mu_1 < mu_2) {
				var truth = "Y";
			} else {
				var truth = "X";
			}
			for (var i = 0;i < samples;i++) {
				var S_x = 0;
				var S_y = 0;
				var S_x2 = 0;
				var S_y2 = 0;
				var n_x = 0;
				var n_y = 0;
				var finished = false;
				while (!finished) {
					// generate samples
					var data_1 = generate([mu_1,var_value]);
					var data_2 = generate([mu_2,var_value]);
					n_x += 1;
					n_y += 1;
					S_x += data_1;
					S_y += data_2;
					S_x2 += data_1*data_1;
					S_y2 += data_2*data_2;
					// test
					var res = checkTest(S_x, S_y, S_x2, S_y2, n_x, n_y)
					if (res) {
						if (res != truth) {
							errs += 1;
						}
						finished = true;
					}
				}
			}
			return errs/samples;
		}*/
	}

	// private functions

	var normal_pac_LR_H0 = function(S_x, S_y, S_x2, S_y2, n_x, n_y, var_value) {
		if (n_x <= 1 || n_y <= 1) {
			return 1;
		}
		var mle_mean = (S_x + S_y)/(n_x + n_y);
		var mle_x = S_x/n_x;
		var mle_y = S_y/n_y;

		var likRatio = Math.exp((n_x+n_y)/(2*2*var_value) * (0.5*mle_x*mle_x + 0.5*mle_y*mle_y - mle_x*mle_y));
		return likRatio;
	}

	functions['normal_pac'] = {
		'l_an' : normal_pac_LR_H0
	}
		  
	tests['normal_pac'] = normal_pac;