	
	/*** test for comparing normal means, known variance ***/

	var normal_known_var_test = function(sides, indifference, type1_error, type2_error, variance) {

		var b0, b1, stoppingTime;

		var x_data = [];
		var y_data = [];
		var n = 0;
		var alpha_value = type1_error;
		var beta_value = type2_error;
		var indiff = indifference;
		var var_value = variance;
		// sufficient stats for test is sum(x_i), sum(y_i), sum(x_i^2) and sum(y_i^2)
		var S_x = 0;
		var S_y = 0;
		var S_x2 = 0;
		var S_y2 = 0;
		var finished = false;
		var L_an;

		/** public functions **/

		this.getResults = function() {
			var L_an = LikH0(S_x, S_y, S_x2, S_y2, n, indiff);
			var L_bn = LikHA(S_x, S_y, S_x2, S_y2, n, indiff);
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
			var std = params[1];
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

			var L_an = LikH0(S_x, S_y, S_x2, S_y2, n, d);
			if (L_an >= b0) {
				return ['false',L_an];
			}
			var L_bn = LikHA(S_x, S_y, S_x2, S_y2, n, d);
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

		var LikH0 = functions['normal_uv'][sides]['l_an'];
		this.LikH0 = LikH0;
		var LikHA = functions['normal_uv'][sides]['l_bn'];
		this.LikHA = LikHA;

		/////////////////////// not checked below

		var boundaryFun = function(indiff) {
			// simulate alpha and beta-value
			var outfun = function(boundaries, n) {
				// calculate alpha with these boundaries
				var results_alpha = alpha(boundaries[0], boundaries[1], indiff, var_bound, simulateResult, n);
				// calculate beta with these boundaries
				var results_beta = beta(boundaries[0], boundaries[1], indiff, var_bound, simulateResult, n);
				return [results_alpha, results_beta];
			}
			return outfun;
		}

		var alpha = functions['normal_uv'][sides]['alpha'];
		this.alpha = alpha;
		var beta = functions['normal_uv'][sides]['beta'];
		this.beta = beta;

		// initialization:
		  // calculate thresholds (unless they are stored in table)
		if (sides in thresholds['normal_uv'] && alpha_value in thresholds['normal_uv'][sides] && beta_value in thresholds['normal_uv'][sides][alpha_value] && indifference in thresholds['normal_uv'][sides][alpha_value][beta_value]) {
			b0 = thresholds['normal_uv'][sides][alpha_value][beta_value][indifference][var_bound][0];
			b1 = thresholds['normal_uv'][sides][alpha_value][beta_value][indifference][var_bound][1];
		} else {
			// calculate thresholds
			console.log("calculating thresholds via simulation")
			//var thr = optimize2d([alpha_value, beta_value], boundaryFun(indifference), [50,10], 0.001, 46000, 400000, 6, 1)
			//var thr = optimize2d([alpha_value, beta_value], boundaryFun(indifference), [98,14.5], 0.001, 46000, 1500000, 6, 1)
			var thr = optimize2d([alpha_value, beta_value], boundaryFun(indifference), [10.9,10.9], 0.001, 46000, 1500000, 6, 1)
			b0 = thr[0];
			b1 = thr[1];
		}

		//this.maxSamplesize = functions['normal_uv'][sides]['max_samplesize'](b0,b1,indiff);
		var simulateH0 = functions['normal_uv'][sides]['simulateH0'](simulateResult, indiff, b0, b1, var_bound);

		// get test variables
		this.properties = {
			'alpha' : alpha_value,
			'beta' : beta_value,
			'indifference region' : indiff,
			'sides' : sides,
			'b0' : b0,
			'b1' : b1,
			'variance bound' : var_bound
		}
	}

	// private functions

	//////////// not checked above

	// change var_bound -> variance
	var normal_kv_twosided_alpha = function(b0, b1, indiff, var_bound, simulateResult, samples) {
		if (!samples) samples = 10000;
		var alphas = [];
		for (var i = 0;i < samples;i++) {
			var res = simulateResult([0,var_bound],[0,var_bound],b0,b1);
			if (res[0] == 'false') {
				alphas.push(1);
			} else {
				alphas.push(0);
			}
		}
		return alphas;
		// TODO : should we include std.dev.?
	}

	// change var_bound -> variance
	var normal_kv_twosided_beta = function(b0, b1, indiff, var_bound, simulateResult, samples) {
		if (!samples) samples = 10000;
		var betas = [];
		for (var i = 0;i < samples;i++) {
			var res = simulateResult([0-indiff/2,var_bound],[0+indiff/2,var_bound],b0,b1);
			if (res[0] == 'true') {
				betas.push(1);
			} else {
				betas.push(0);
			}
		}
		return betas;
		// TODO : should we include std.dev.?
	}

	//changed
	var normal_kv_twosided_LR_H0 = function(S_x, S_y, S_x2, S_y2, n, indiff) {
		if (n == 1) {
			return 1;
		}
		var mle_mean = (S_x + S_y)/(2*n);
		var mle_x = S_x/n;
		var mle_y = S_y/n;

		var likRatio = Math.exp(n/(2*var_value) * (0.5*mle_x*mle_x + 0.5*mle_y*mle_y - mle_x*mle_y));
		return likRatio;
	}

	// changed
	var normal_kv_twosided_LR_HA = function(S_x, S_y, S_x2, S_y2, n, indiff) {
		if (n == 1) {
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

		var mle_lik_part = S_x2 - n*unc_mle_x*unc_mle_x + S_y2 - n*unc_mle_y*unc_mle_y;
		if (pos_lik_part < neg_lik_part) {
			return Math.exp( 1/(2*var_value)*(pos_lik_part - mle_lik_part) );
		} else {
			return Math.exp( 1/(2*var_value)*(neg_lik_part - mle_lik_part) );
		}
	}

	// change var_bound -> variance
	var normal_uv_twosided_simulateH0 = function(simRes, indiff, b0, b1, var_bound) {
		var returnFun = function() {
			var res = simRes([0,var_bound],[0,var_bound],b0,b1)[4];
			return res;
		}
		return returnFun;
	}

	// change var_bound -> variance
	var normal_kv_onesided_alpha = function(b0, b1, indiff, var_bound, simulateResult, samples) {
		if (!samples) samples = 10000;
		var alphas = [];
		for (var i = 0;i < samples;i++) {
			var res = simulateResult([0-indiff/2,var_bound],[0+indiff/2,var_bound],b0,b1);
			if (res[0] == 'false') {
				alphas.push(1);
			} else {
				alphas.push(0);
			}
		}
		return alphas;
	}

	// change var_bound -> variance
	var normal_kv_onesided_beta = function(b0, b1, indiff, var_bound, simulateResult, samples) {
		if (!samples) samples = 10000;
		var betas = [];
		for (var i = 0;i < samples;i++) {
			var res = simulateResult([0+indiff/2,var_bound],[0-indiff/2,var_bound],b0,b1);
			if (res[0] == 'true') {
				betas.push(1);
			} else {
				betas.push(0);
			}
		}
		return betas;
		// TODO : should we include std.dev.?
	}

	// changed
	var normal_kv_onesided_LR_H0 = function(S_x, S_y, S_x2, S_y2, n, indiff) {
		// H0 is that mu_1 < mu_2 -indiff -> mu_1 = mu_2 - indiff
		if (n == 1) {
			return 1;
		}
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

	// changed
	var normal_uv_onesided_LR_HA = function(S_x, S_y, S_x2, S_y2, n, indiff) {
		if (n == 1) {
			return 1;
		}
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

	// change car_bound -> variance
	var normal_uv_onesided_simulateH0 = function(simRes, indiff, b0, b1, var_bound) {
		var returnFun = function() {
			var res = simRes([0-indiff/2,var_bound],[0+indiff/2,var_bound],b0,b1)[4];
			return res;
		}
		return returnFun;
	}
	//////////// not checked below

	functions['normal_uv'] = {
		'two-sided' : {
			'l_an' : normal_uv_twosided_LR_H0,
			'l_bn' : normal_uv_twosided_LR_HA,
			'alpha' : normal_uv_twosided_alpha,
			'beta' : normal_uv_twosided_beta,
			'simulateH0' : normal_uv_twosided_simulateH0,
		},
		'one-sided' : {
			'l_an' : normal_uv_onesided_LR_H0,
			'l_bn' : normal_uv_onesided_LR_HA,
			'alpha' : normal_uv_onesided_alpha,
			'beta' : normal_uv_onesided_beta,
			'simulateH0' : normal_uv_onesided_simulateH0,	
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
						1 : [135, 135], // variance bound
					}
				}
			}
		}
	}

	tests['normal_unknown_var'] = normal_unknown_var_test;