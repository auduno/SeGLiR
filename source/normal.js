	
	/*** test for comparing normal means, unknown or known variance ***/

	var normal_test = function(sides, indifference, type1_error, type2_error, variance, variance_bound) {

		var b0, b1, stoppingTime, var_bound, var_value;

		// check input
		if (sides != "one-sided" && sides != "two-sided") {
			console.log("parameter 'sides' must be either 'one-sided' or 'two-sided', input was : '"+sides+"'!");
			return;
		}
		if (typeof(indifference) != 'number' || indifference <= 0) {
			console.log("parameter 'indifference' must be a number above zero, input was : "+indifference);
			return;
		}
		if (typeof(type1_error) != 'number' || type1_error <= 0 || type1_error >= 1) {
			console.log("parameter 'type1_error' must be a number between 0 and 1, input was : "+type1_error);
			return;
		}
		if (typeof(type2_error) != 'number' || type2_error <= 0 || type2_error >= 1) {
			console.log("parameter 'type2_error' must be a number between 0 and 1, input was : "+type2_error);
			return;
		}
		if (typeof(variance) == 'undefined') {
			if (typeof(variance_bound) != 'number' || variance_bound <= 0) {
				console.log("when parameter 'variance' is undefined, 'variance_bound' must be a valid variance, i.e. number above 0, input was : "+variance_bound);
				return;
			}
		} else if (typeof(variance) != 'number' || variance <= 0) {
			console.log("when parameter 'variance' is specified, it must be a valid variance, i.e. number above 0, input was : "+variance);
			return;
		}

		var x_data = [];
		var y_data = [];
		var n = 0;
		var alpha_value = type1_error;
		var beta_value = type2_error;
		var indiff = indifference;
		var var_bound 
		if (typeof(variance) == "undefined") {
			var_bound = variance_bound;
		} else {
			var_value = variance;
		}
		// sufficient stats for test is sum(x_i), sum(y_i), sum(x_i^2) and sum(y_i^2)
		var S_x = 0;
		var S_y = 0;
		var S_x2 = 0;
		var S_y2 = 0;
		var finished = false;
		var L_an;

		/** public functions **/

		this.getResults = function() {
			var L_an = LikH0(S_x, S_y, S_x2, S_y2, n, indiff, var_value);
			var L_bn = LikHA(S_x, S_y, S_x2, S_y2, n, indiff, var_value);
			return {
				'S_x' : S_x,
				'S_y' : S_y,
				'S_x2' : S_x2,
				'S_y2' : S_y2,
				'L_an' : L_an,
				'L_bn' : L_bn,
				'finished' : finished,
				'n' : n
			};
		}

		// get p-value (only when test is done)
		this.pValue = function(samples) {
			if (!finished) {
				return undefined;
			}
			if (!samples) samples = 10000;
			console.log("calculating p-value via simulation");
			var res = 0;
			for (var i = 0;i < samples;i++) {
				if (simulateH0() >= L_an) {
					res += 1;
				}
			}
			return res/samples;
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
				var res = simulateResult(ests[0],ests[1],b0,b1);
				var time = res[5];
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
		this.addData = function(points) {
			if (finished) {
				if (typeof points[0] === 'number') x_data.push(points[0]);
				if (typeof points[1] === 'number') y_data.push(points[1]);
			} else {
				if (typeof points[0] === 'number' && typeof points[1] === 'number') {
					if (x_data.length == y_data.length) {
						S_x += points[0];
						S_y += points[1];
						S_x2 += points[0]*points[0];
						S_y2 += points[1]*points[1];
						n += 1;
					} else if (x_data.length > y_data.length) {
						S_x += x_data[n];
						S_y += points[1];
						S_x2 += x_data[n]*x_data[n];
						S_y2 += points[1]*points[1];
						n += 1;
					} else {
						S_x += points[0];
						S_y += y_data[n];
						S_x2 += points[0]*points[0];
						S_y2 += y_data[n]*y_data[n];
						n += 1;
					}
					x_data.push(points[0])
					y_data.push(points[1])
				} else if (typeof points[0] === 'number') {
					if (x_data.length < y_data.length) {
						S_x += points[0];
						S_y += y_data[n];
						S_x2 += points[0]*points[0];
						S_y2 += y_data[n]*y_data[n];
						n += 1;
					}
					x_data.push(points[0]);
				} else if (typeof points[1] === 'number') {
					if (x_data.length > y_data.length) {
						S_x += x_data[n];
						S_y += points[1];
						S_x2 += x_data[n]*x_data[n];
						S_y2 += points[1]*points[1];
						n += 1;
					}
					y_data.push(points[1]);
				}
			}
			
			var result = checkTest(S_x, S_y, S_x2, S_y2, n, indiff, b0, b1);
			if (result) {
				finished = true;
				stoppingTime = n;
				L_an = result[1];
				return result[0];
			}
		}

		// get expected samplesize for some parameters
		this.expectedSamplesize = function(params_1, params_2, samples) {
			// simulate it enough times
			if (!samples) samples = 10000;
			console.log("calculating expected samplesize via simulation");
			var times = [];
			for (var i = 0;i < samples;i++) {
				var res = simulateResult(params_1,params_2,b0,b1)
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

		var simulateResult = function(params_1, params_2, b0, b1) {
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
				var result = checkTest(S_x, S_y, S_x2, S_y2, time, indiff, b0, b1);
				if (result) finished = true;
			}
			// return result, S_x, S_y, stoppingTime
			return [result[0], S_x, S_y, S_x2, S_y2, time, result[1]];
		}
		this.simulateResult = simulateResult;

		var checkTest = function(S_x, S_y, S_x2, S_y2, n, d, b0, b1) {
			// check if test should be stopped
			
			// TODO : should I check for when both L_an and L_bn pass thresholds?

			var L_an = LikH0(S_x, S_y, S_x2, S_y2, n, d, var_value);
			if (L_an >= b0) {
				return ['false',L_an];
			}
			var L_bn = LikHA(S_x, S_y, S_x2, S_y2, n, d, var_value);
			if (L_bn >= b1) {
				return ['true',L_an]
			}
			return undefined
		}
		
		var biasFun = function() {
			var outfun = function(parameters, n) {
				var results_p1 = []
				var results_p2 = []
				for (var i = 0;i < n;i++) {
					// generate sequences
					var res = simulateResult(parameters[0], parameters[1], b0, b1);
					results_p1.push( res[1]/res[5] );
					results_p2.push( res[2]/res[5] );
				}
				return [results_p1, results_p2];
			}
			return outfun;
		}

		if (var_value) {
			var LikH0 = functions['normal_kv'][sides]['l_an'];
			var LikHA = functions['normal_kv'][sides]['l_bn'];
		} else {
			var LikH0 = functions['normal_uv'][sides]['l_an'];
			var LikHA = functions['normal_uv'][sides]['l_bn'];
		}
		this.LikH0 = LikH0;
		this.LikHA = LikHA;

		var boundaryFun = function(indiff) {
			// simulate alpha and beta-value
			var outfun = function(boundaries, n) {
				if (var_value) {
					var results_alpha = alpha(boundaries[0], boundaries[1], indiff, var_value, simulateResult, n);
					var results_beta = beta(boundaries[0], boundaries[1], indiff, var_value, simulateResult, n);
				} else {
					var results_alpha = alpha(boundaries[0], boundaries[1], indiff, var_bound, simulateResult, n);
					var results_beta = beta(boundaries[0], boundaries[1], indiff, var_bound, simulateResult, n);
				}
				return [results_alpha, results_beta];
			}
			return outfun;
		}

		if (var_value) {
			var alpha = functions['normal_kv'][sides]['alpha'];
			var beta = functions['normal_kv'][sides]['beta'];
		} else {
			var alpha = functions['normal_uv'][sides]['alpha'];
			var beta = functions['normal_uv'][sides]['beta'];
		}
		this.alpha = alpha;
		this.beta = beta;

		// initialization:
		  // calculate thresholds (unless they are stored in table)

		if (var_value) {
			var our_thresholds = thresholds['normal_kv'];
			var our_var = var_value;
		} else {
			var our_thresholds = thresholds['normal_uv'];
			var our_var = var_bound;
		}
		if (sides in our_thresholds && alpha_value in our_thresholds[sides] && beta_value in our_thresholds[sides][alpha_value] && indifference in our_thresholds[sides][alpha_value][beta_value] && our_var in our_thresholds[sides][alpha_value][beta_value][indifference]) {
			b0 = our_thresholds[sides][alpha_value][beta_value][indifference][our_var][0];
			b1 = our_thresholds[sides][alpha_value][beta_value][indifference][our_var][1];
		} else {
			// calculate thresholds
			console.log("calculating thresholds via simulation")
			//var thr = optimize2d([alpha_value, beta_value], boundaryFun(indifference), [50,10], 0.001, 46000, 400000, 6, 1)
			//var thr = optimize2d([alpha_value, beta_value], boundaryFun(indifference), [98,14.5], 0.001, 46000, 1500000, 6, 1)
			var thr = optimize2d([alpha_value, beta_value], boundaryFun(indifference), [200,9], 0.001, 46000, 1500000, 6, 1)
			b0 = thr[0];
			b1 = thr[1];
		}

		// TODO : implement this for known variance
		//this.maxSamplesize = functions['normal_uv'][sides]['max_samplesize'](b0,b1,indiff);
		if (var_value) {
			var simulateH0 = functions['normal_kv'][sides]['simulateH0'](simulateResult, indiff, b0, b1, var_value);
		} else {
			var simulateH0 = functions['normal_uv'][sides]['simulateH0'](simulateResult, indiff, b0, b1, var_bound);
		}

		// get test variables
		this.properties = {
			'alpha' : alpha_value,
			'beta' : beta_value,
			'indifference region' : indiff,
			'sides' : sides,
			'b0' : b0,
			'b1' : b1,
			'variance' : var_value,
			'variance bound' : var_bound
		}
	}

	// private functions

	var normal_twosided_alpha = function(b0, b1, indiff, var_val, simulateResult, samples) {
		if (!samples) samples = 10000;
		var alphas = [];
		//var starttime = (new Date()).getTime();
		for (var i = 0;i < samples;i++) {
			var res = simulateResult([0,var_val],[0,var_val],b0,b1);
			if (res[0] == 'false') {
				alphas.push(1);
			} else {
				alphas.push(0);
			}
		}
		//console.log("time:"+( (new Date()).getTime()-starttime ))
		//console.log("mean:"+mean(alphas));
		//console.log("std_err:"+boot_std(alphas,1000));
		return alphas;
	}

	var normal_twosided_alpha_imp = function(b0, b1, indiff, var_val, simulateResult, samples) {
		if (!samples) samples = 10000;
		var alphas = [];
		var beta = 1; // precision/inverse-variance of the importance sampling distribution
		//var starttime = (new Date()).getTime();
		for (var i = 0;i < samples;i++) {
			var z = jStat.jStat.normal.sample(0,Math.sqrt(1/beta));

			var finished = false;
			var S_x = 0;
			var n = 0;
			var result = undefined;
			while (!finished) {
				n += 1;
				// pull xs from N(0,2*var_val)
				S_x += jStat.jStat.normal.sample(z,Math.sqrt(2*var_val));
				// test on simplified boundaries
				var L_na = Math.exp( S_x*S_x/(4*n*var_val) );
				if (L_na >= b0) {
					finished = true;
					result = "false"
				}
				if (Math.abs(S_x/n) < indiff) {
					if (S_x/n > 0) {
						var L_nb = Math.exp( n*(S_x/n - indiff)*(S_x/n - indiff)/(4*var_val) );
					} else {
						var L_nb = Math.exp( n*(S_x/n + indiff)*(S_x/n + indiff)/(4*var_val) );
					}
					if (L_nb >= b1) {
						finished = true;
						result = "true";
					}
				}
			}

			if (result == 'false') {
				var b2v = 2*beta*var_val;
				var weight = Math.sqrt(b2v/(n + b2v))*Math.exp( (S_x*S_x)/(4*var_val*(n+b2v)) );
				alphas.push(1/weight);
			} else {
				alphas.push(0);
			}
		}
		//console.log("time:"+( (new Date()).getTime()-starttime ))
		//console.log("mean:"+mean(alphas));
		//console.log("std_err:"+boot_std(alphas,1000));
		return alphas;
	}

	var normal_twosided_beta = function(b0, b1, indiff, var_val, simulateResult, samples) {
		if (!samples) samples = 10000;
		var betas = [];
		//var starttime = (new Date()).getTime();
		for (var i = 0;i < samples;i++) {
			var res = simulateResult([-indiff/2,var_val],[indiff/2,var_val],b0,b1);
			if (res[0] == 'true') {
				betas.push(1);
			} else {
				betas.push(0);
			}
		}
		//console.log("time:"+( (new Date()).getTime()-starttime ))
		//console.log("mean:"+mean(betas));
		//console.log("std_err:"+boot_std(betas,1000));
		return betas;
	}

	var normal_twosided_beta_imp = function(b0, b1, indiff, var_val, simulateResult, samples) {
		if (!samples) samples = 10000;
		var betas = [];
		//var starttime = (new Date()).getTime();
		for (var i = 0;i < samples;i++) {
			var finished = false;
			var S_x = 0;
			var n = 0;
			var result = undefined;
			while (!finished) {
				n += 1;
				// pull xs from N(0,2*var_val)
				S_x += jStat.jStat.normal.sample(0,Math.sqrt(2*var_val));
				// test on simplified boundaries
				var L_na = Math.exp( S_x*S_x/(4*n*var_val) );
				if (L_na >= b0) {
					finished = true;
					result = "false"
				}
				if (Math.abs(S_x/n) < indiff) {
					if (S_x/n > 0) {
						var L_nb = Math.exp( n*(S_x/n - indiff)*(S_x/n - indiff)/(4*var_val) );
					} else {
						var L_nb = Math.exp( n*(S_x/n + indiff)*(S_x/n + indiff)/(4*var_val) );
					}
					//var L_nb = Math.exp( (S_x + n*indiff)*(S_x + n*indiff)/(4*n*var_val) );
					if (L_nb >= b1) {
						finished = true;
						result = "true";
					}
				}
			}

			if (result == 'true') {
				var weight = Math.exp( (2*indiff*S_x - n*indiff*indiff)/(4*var_val) );
				betas.push(weight);
			} else {
				betas.push(0);
			}
		}
		//console.log("time:"+( (new Date()).getTime()-starttime ))
		//console.log("mean:"+mean(betas));
		//console.log("std_err:"+boot_std(betas,1000));
		return betas;
	}

	var normal_uv_twosided_LR_H0 = function(S_x, S_y, S_x2, S_y2, n, indiff) {
		if (n == 1) {
			// TODO : when n = 1, mle_llik is always -Infinity, how to deal with this?
			return 1;
		}
		var mle_mean = (S_x + S_y)/(2*n);
		var part1 = S_x2 + S_y2 - 2*n*mle_mean*mle_mean
		if (part1 == 0) {
			var likRatio = 0;
		} else {
			var likRatio = Math.exp(n * ( Math.log( part1 ) - Math.log( S_x2 - n*(S_x/n)*(S_x/n) + S_y2 - n*(S_y/n)*(S_y/n) ) ));
		}

		return likRatio;
	}

	var normal_uv_twosided_LR_HA = function(S_x, S_y, S_x2, S_y2, n, indiff) {
		if (n == 1) {
			// TODO : when n = 1, mle_llik is always -Infinity, how to deal with this?
			return 1;
		}

		var unc_mle_x = S_x/n;
		var unc_mle_y = S_y/n;

		if (Math.abs(unc_mle_x-unc_mle_y) > indiff) {
			return 1;
		}
		
		var pos = 0.5*(unc_mle_x + unc_mle_y + indiff); // mle of mu_1 when constrained so mu_1 = mu_2 + d
		var neg = 0.5*(unc_mle_x + unc_mle_y - indiff); // mle of mu_1 when constrained so mu_1 = mu_2 - d

		var pos_lik_part = S_x2 - 2*pos*S_x + n*pos*pos + S_y2 - 2*S_y*(pos-indiff) + n*(pos-indiff)*(pos-indiff);
		var neg_lik_part = S_x2 - 2*neg*S_x + n*neg*neg + S_y2 - 2*S_y*(neg+indiff) + n*(neg+indiff)*(neg+indiff);

		var mle_lik_part = Math.log(S_x2 - n*unc_mle_x*unc_mle_x + S_y2 - n*unc_mle_y*unc_mle_y );
		if (pos_lik_part < neg_lik_part) {
			return Math.exp( n*(Math.log(pos_lik_part) - mle_lik_part) );
		} else {
			return Math.exp( n*(Math.log(neg_lik_part) - mle_lik_part) );
		}
	}

	var normal_twosided_simulateH0 = function(simRes, indiff, b0, b1, var_val) {
		var returnFun = function() {
			var res = simRes([0,var_val],[0,var_val],b0,b1)[4];
			return res;
		}
		return returnFun;
	}

	var normal_onesided_alpha = function(b0, b1, indiff, var_val, simulateResult, samples) {
		if (!samples) samples = 10000;
		var alphas = [];
		for (var i = 0;i < samples;i++) {
			var res = simulateResult([0-indiff/2,var_val],[0+indiff/2,var_val],b0,b1);
			if (res[0] == 'false') {
				alphas.push(1);
			} else {
				alphas.push(0);
			}
		}
		return alphas;
	}

	var normal_onesided_beta = function(b0, b1, indiff, var_val, simulateResult, samples) {
		if (!samples) samples = 10000;
		var betas = [];
		for (var i = 0;i < samples;i++) {
			var res = simulateResult([0+indiff/2,var_val],[0-indiff/2,var_val],b0,b1);
			if (res[0] == 'true') {
				betas.push(1);
			} else {
				betas.push(0);
			}
		}
		return betas;
		// TODO : should we include std.dev.?
	}

	var normal_uv_onesided_LR_H0 = function(S_x, S_y, S_x2, S_y2, n, indiff) {
		if (n == 1) {
			// TODO : when n = 1, mle_llik is always -Infinity, how to deal with this?
			return 1;
		}
		// H0 is that mu_1 < mu_2 -> mu_1 - mu_2 < -indiff -> mu_1 = mu_2 - indiff
		var unc_mle_x = S_x/n;
		var unc_mle_y = S_y/n;

		if (unc_mle_x-unc_mle_y <= -indiff) {
			return 1;
		}
		
		var neg = 0.5*(unc_mle_x + unc_mle_y - indiff); // mle of mu_1 when constrained so mu_1 = mu_2 - d -> mu_1 <= mu_2 - d -> mu_1 - mu_2 <= - d
		var neg_llik = S_x2 - 2*neg*n*unc_mle_x + n*neg*neg + S_y2 - 2*n*unc_mle_y*(neg+indiff) + n*(neg+indiff)*(neg+indiff);
		var mle_llik = Math.log(S_x2 - n*unc_mle_x*unc_mle_x + S_y2 - n*unc_mle_y*unc_mle_y );
		
		return Math.exp( n*(Math.log(neg_llik) - mle_llik) );
	}

	var normal_uv_onesided_LR_HA = function(S_x, S_y, S_x2, S_y2, n, indiff) {
		if (n == 1) {
			// TODO : when n = 1, mle_llik is always -Infinity, how to deal with this?
			return 1;
		}
		var unc_mle_x = S_x/n;
		var unc_mle_y = S_y/n;

		if (unc_mle_x-unc_mle_y >= indiff) {
			return 1;
		}
		
		var pos = 0.5*(unc_mle_x + unc_mle_y + indiff); // mle of mu_1 when constrained so mu_1 = mu_2 + d -> mu_1 >= mu_2 + d -> mu_1 - mu_2 >= d
		var pos_llik = S_x2 - 2*pos*n*unc_mle_x + n*pos*pos + S_y2 - 2*n*unc_mle_y*(pos-indiff) + n*(pos-indiff)*(pos-indiff);
		var mle_llik = Math.log(S_x2 - n*unc_mle_x*unc_mle_x + S_y2 - n*unc_mle_y*unc_mle_y );

		return Math.exp( n*(Math.log(pos_llik) - mle_llik) );
	}

	var normal_onesided_simulateH0 = function(simRes, indiff, b0, b1, var_val) {
		var returnFun = function() {
			var res = simRes([0-indiff/2,var_val],[0+indiff/2,var_val],b0,b1)[4];
			return res;
		}
		return returnFun;
	}

	var normal_kv_twosided_LR_H0 = function(S_x, S_y, S_x2, S_y2, n, indiff, var_value) {
		var mle_mean = (S_x + S_y)/(2*n);
		var mle_x = S_x/n;
		var mle_y = S_y/n;

		var likRatio = Math.exp(n/(2*var_value) * (0.5*mle_x*mle_x + 0.5*mle_y*mle_y - mle_x*mle_y));
		//var likRatio2 = Math.exp((S_x-S_y)*(S_x-S_y)/(4*n*var_value));

		return likRatio;
	}

	var normal_kv_twosided_LR_HA = function(S_x, S_y, S_x2, S_y2, n, indiff, var_value) {
		var unc_mle_x = S_x/n;
		var unc_mle_y = S_y/n;

		if (Math.abs(unc_mle_x-unc_mle_y) > indiff) {
			return 1;
		}
		
		var pos = 0.5*(unc_mle_x + unc_mle_y + indiff); // mle of mu_1 when constrained so mu_1 = mu_2 + d
		var neg = 0.5*(unc_mle_x + unc_mle_y - indiff); // mle of mu_1 when constrained so mu_1 = mu_2 - d

		var pos_lik_part = S_x2 - 2*pos*S_x + n*pos*pos + S_y2 - 2*S_y*(pos-indiff) + n*(pos-indiff)*(pos-indiff);
		var neg_lik_part = S_x2 - 2*neg*S_x + n*neg*neg + S_y2 - 2*S_y*(neg+indiff) + n*(neg+indiff)*(neg+indiff);

		var mle_lik_part = S_x2 - n*unc_mle_x*unc_mle_x + S_y2 - n*unc_mle_y*unc_mle_y;
		if (pos_lik_part < neg_lik_part) {
			//console.log("regular pos likratio : "+Math.exp( 1/(2*var_value)*(pos_lik_part - mle_lik_part) ) );
			//console.log("simple pos likratio : "+Math.exp( (S_x-S_y-n*indiff)*(S_x-S_y-n*indiff)/(4*n*var_value) ) );
			return Math.exp( 1/(2*var_value)*(pos_lik_part - mle_lik_part) );
		} else {
			//console.log("regular neg likratio : "+Math.exp( 1/(2*var_value)*(neg_lik_part - mle_lik_part) ) );
			//console.log("simple neg likratio : "+Math.exp( (S_x-S_y+n*indiff)*(S_x-S_y+n*indiff)/(4*n*var_value) ) );
			return Math.exp( 1/(2*var_value)*(neg_lik_part - mle_lik_part) );
		}
	}

	var normal_kv_onesided_LR_H0 = function(S_x, S_y, S_x2, S_y2, n, indiff, var_value) {
		// H0 is that mu_1 < mu_2 -indiff -> mu_1 = mu_2 - indiff
		var unc_mle_x = S_x/n;
		var unc_mle_y = S_y/n;

		if (unc_mle_x-unc_mle_y <= -indiff) {
			return 1;
		}
		
		var neg = 0.5*(unc_mle_x + unc_mle_y - indiff); // mle of mu_1 when constrained so mu_1 = mu_2 - d -> mu_1 <= mu_2 - d -> mu_1 - mu_2 <= - d
		var neg_lik_part = S_x2 - 2*neg*S_x + n*neg*neg + S_y2 - 2*S_y*(neg+indiff) + n*(neg+indiff)*(neg+indiff);
		var mle_lik_part = S_x2 - n*unc_mle_x*unc_mle_x + S_y2 - n*unc_mle_y*unc_mle_y;

		return Math.exp( 1/(2*var_value)*(neg_lik_part - mle_lik_part) );
	}

	var normal_kv_onesided_LR_HA = function(S_x, S_y, S_x2, S_y2, n, indiff, var_value) {
		var unc_mle_x = S_x/n;
		var unc_mle_y = S_y/n;

		if (unc_mle_x-unc_mle_y >= indiff) {
			return 1;
		}
		
		var pos = 0.5*(unc_mle_x + unc_mle_y + indiff); // mle of mu_1 when constrained so mu_1 = mu_2 + d -> mu_1 >= mu_2 + d -> mu_1 - mu_2 >= d
		var pos_lik_part = S_x2 - 2*pos*S_x + n*pos*pos + S_y2 - 2*S_y*(pos-indiff) + n*(pos-indiff)*(pos-indiff);
		var mle_lik_part = S_x2 - n*unc_mle_x*unc_mle_x + S_y2 - n*unc_mle_y*unc_mle_y;

		return Math.exp( 1/(2*var_value)*(pos_lik_part - mle_lik_part) );
	}

	functions['normal_uv'] = {
		'two-sided' : {
			'l_an' : normal_uv_twosided_LR_H0,
			'l_bn' : normal_uv_twosided_LR_HA,
			'alpha' : normal_twosided_alpha,
			'beta' : normal_twosided_beta,
			'simulateH0' : normal_twosided_simulateH0,
		},
		'one-sided' : {
			'l_an' : normal_uv_onesided_LR_H0,
			'l_bn' : normal_uv_onesided_LR_HA,
			'alpha' : normal_onesided_alpha,
			'beta' : normal_onesided_beta,
			'simulateH0' : normal_onesided_simulateH0,	
		}
	}
	functions['normal_kv'] = {
		'two-sided' : {
			'l_an' : normal_kv_twosided_LR_H0,
			'l_bn' : normal_kv_twosided_LR_HA,
			'alpha' : normal_twosided_alpha_imp,
			'beta' : normal_twosided_beta_imp,
			'simulateH0' : normal_twosided_simulateH0,
		},
		'one-sided' : {
			'l_an' : normal_kv_onesided_LR_H0,
			'l_bn' : normal_kv_onesided_LR_HA,
			'alpha' : normal_onesided_alpha,
			'beta' : normal_onesided_beta,
			'simulateH0' : normal_onesided_simulateH0,	
		}
	}
		  
	// precalculated thresholds
	thresholds['normal_uv'] = {
		'two-sided' : {
			0.05 : { // alpha
				0.10 : { // beta
					0.1 : { // indifference
						1 : [430, 8.85], // variance bound
					}
				}
			}
		},
		'one-sided' : {
			0.05 : { // alpha
				0.05 : { // beta
					0.1 : { // indifference
						1 : [135, 135], // check result
					}
				}
			}
		}
	}
	thresholds['normal_kv'] = {
		'two-sided' : {
			0.05 : { // alpha
				0.10 : { // beta
					0.2 : { // indifference
						1 : [106, 8.0] // variance
					},
					0.1 : { // indifference
						1 : [137.5, 8.3] // variance
					},
					0.05 : { // indifference
						1 : [170.5, 8.55] // variance
					},
					/*0.025 : { // indifference
						1 : [200, 9] // check result
					}*/
				}
			}
		},
		'one-sided' : {
			0.05 : { // alpha
				0.05 : { // beta
					0.1 : { // indifference
						1 : [49, 49], // variance, TODO : check result
					}
				}
			}
		}
	}

	tests['normal'] = normal_test;