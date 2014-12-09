	
	/*** test for bernoulli proportions ***/

	var bernoulli_test = function(sides, indifference, type1_error, type2_error, simulateThreshold) {

		var b0, b1, stoppingTime;

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
		if (typeof(simulateThreshold) == "undefined") {
			simulateThreshold = true;
		}

		var x_data = [];
		var y_data = [];
		var n = 0;
		var alpha_value = type1_error;
		var beta_value = type2_error;
		var indiff = indifference;
		var S_x = 0;
		var S_y = 0;
		var finished = false;
		var L_an;

		/** public functions **/

		this.getResults = function() {
			var L_an = LikH0(S_x, S_y, n, indiff);
			var L_bn = LikHA(S_x, S_y, n, indiff);
			return {
				'S_x' : S_x,
				'S_y' : S_y,
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
		this.estimate = function(max_samples) {
			if (!finished) {
				return undefined;
			}
			if (typeof(max_samples) == "undefined") {
				max_samples = 1500000;
			}
			var ests = optimize2d([S_x/n, S_y/n], biasFun(), [S_x/n, S_y/n], 0.005, 16400, max_samples, 0.02, 0, 1, true);
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
			if (!simulateThreshold) {
				console.log("No thresholds are defined, this mode is only for manually finding thresholds.")
				return;
			}
			if (finished) {
				if (typeof points['x'] === 'number') x_data.push(points['x']);
				if (typeof points['y'] === 'number') y_data.push(points['y']);
			} else {
				if (typeof points['x'] === 'number' && typeof points['y'] === 'number') {
					if (x_data.length == y_data.length) {
						S_x += points['x'];
						S_y += points['y'];
						n += 1;
					} else if (x_data.length > y_data.length) {
						S_y += points['y'];
						S_x += x_data[n];
						n += 1;
					} else {
						S_x += points['x'];
						S_y += y_data[n];
						n += 1;
					}
					x_data.push(points['x'])
					y_data.push(points['y'])
				} else if (typeof points['x'] === 'number') {
					if (x_data.length < y_data.length) {
						S_x += points['x'];
						S_y += y_data[n];
						n += 1;
					}
					x_data.push(points['x']);
				} else if (typeof points['y'] === 'number') {
					if (x_data.length > y_data.length) {
						S_y += points['y'];
						S_x += x_data[n];
						n += 1;
					}
					y_data.push(points['y']);
				}
			}
			
			var result = checkTest(S_x, S_y, n, indiff, b0, b1);
			if (result) {
				finished = true;
				stoppingTime = n;
				L_an = result[1];
				return result[0];
			}
		}

		// get expected samplesize for some parameters
		this.expectedSamplesize = function(p1, p2, samples) {
			if (!simulateThreshold) {
				console.log("No thresholds are defined, this mode is only for manually finding thresholds.")
				return;
			}
			// simulate it enough times
			if (!samples) samples = 10000;
			console.log("calculating expected samplesize via simulation");
			var times = [];
			for (var i = 0;i < samples;i++) {
				var res = simulateResult(p1,p2,b0,b1)
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
					var res = simulateResult(pt[0], pt[1], b0, b1);
					results_p1.push( res[1]/res[3] );
					results_p2.push( res[2]/res[3] );
				}
				return [results_p1, results_p2];
			}
			return outfun;
		}
		
		var checkTest = function(S_x, S_y, n, d, b0, b1) {
			// check if test should be stopped
			
			// TODO : should I check for cases when both L_an and L_bn pass thresholds?

			var L_an = LikH0(S_x, S_y, n, d);
			if (L_an >= b0) {
				return ['false',L_an];
			}
			var L_bn = LikHA(S_x, S_y, n, d);
			if (L_bn >= b1) {
				return ['true',L_an]
			}
			return undefined
		}

		var LikH0 = functions['bernoulli'][sides]['l_an'];
		var LikHA = functions['bernoulli'][sides]['l_bn'];

		var boundaryFun = function(indiff) {
			// simulate alpha and beta-value
			var outfun = function(boundaries, n) {
				// calculate alpha with these boundaries
				var results_alpha = alpha(boundaries[0], boundaries[1], indiff, simulateResult, n);
				// calculate beta with these boundaries
				var results_beta = beta(boundaries[0], boundaries[1], indiff, simulateResult, n);
				return [results_alpha, results_beta];
			}
			return outfun;
		}

		var generate = function(p) {
			if (Math.random() < p) {return 1;} else {return 0;}
		}
		
		var alpha = functions['bernoulli'][sides]['alpha'];
		var beta = functions['bernoulli'][sides]['beta'];
				
		var simulateResult = function(p1, p2, b0, b1) {
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
				var result = checkTest(S_x, S_y, time, indiff, b0, b1);
				if (result) finished = true;
			}
			// return result, S_x, S_y, stoppingTime
			return [result[0], S_x, S_y, time, result[1]];
		}

		this.alpha_level = function(b0,b1,samples) {
			var alphas = alpha(b0, b1, indiff, simulateResult, samples);
			var mn = mean(alphas);
			var sderr = boot_std(alphas,1000)
			return [mn,sderr];
		}

		this.beta_level = function(b0,b1,samples) {
			var betas = beta(b0, b1, indiff, simulateResult, samples);
			var mn = mean(betas);
			var sderr = boot_std(betas,1000)
			return [mn,sderr];	
		}

		// initialization:
		  // calculate thresholds (unless they are stored in table)
		if (sides in thresholds['bernoulli'] && alpha_value in thresholds['bernoulli'][sides] && beta_value in thresholds['bernoulli'][sides][alpha_value] && indifference in thresholds['bernoulli'][sides][alpha_value][beta_value]) {
			b0 = thresholds['bernoulli'][sides][alpha_value][beta_value][indifference][0];
			b1 = thresholds['bernoulli'][sides][alpha_value][beta_value][indifference][1];
		} else if (simulateThreshold) {
			// calculate thresholds
			console.log("Calculating thresholds via simulation.")
			console.log("Please note : Calculating thresholds via simulation might take a long time. To save time, consult the SeGLiR reference to find test settings that already have precalculated thresholds.")
			//var thr = optimize2d([alpha_value, beta_value], boundaryFun(indifference), [50,10], 0.001, 46000, 400000, 6, 1)
			//var thr = optimize2d([alpha_value, beta_value], boundaryFun(indifference), [98,14.5], 0.001, 46000, 1500000, 6, 1)
			var thr = optimize2d([alpha_value, beta_value], boundaryFun(indifference), [10,10], 0.001, 46000, 1500000, 6, 1, undefined, true, false);
			b0 = thr[0];
			b1 = thr[1];
		} else {
			console.log("NB! No precalculated thresholds are found and simulation of thresholds is disabled - this mode is only for manually finding thresholds for a given alpha- and beta-level.")
		}

		this.maxSamplesize = functions['bernoulli'][sides]['max_samplesize'](b0,b1,indiff);
		var simulateH0 = functions['bernoulli'][sides]['simulateH0'](simulateResult, indiff, b0, b1);

		// get test variables
		this.properties = {
			'alpha' : alpha_value,
			'beta' : beta_value,
			'indifference region' : indiff,
			'sides' : sides,
			'b0' : b0,
			'b1' : b1
		}
	}

	// private functions

	var solveConstrainedBinomialMLE = function(S_x, S_y, n, d) {
		// solves MLE of p1 with the constraint that p1 = p2 - d
		var a = (3*d*n - S_x - S_y - 2*n);
		var b = (S_x - 2*d*S_x + S_y - 2*d*n + d*d*n);
		var P = -a/(6*n);
		var Q = P*P*P + (a*b - 3*2*n*(d*S_x - d*d*S_x))/(6*2*2*n*n);
		var R = b/(6*n);
		var innerSquare = Q*Q + (R - P*P)*(R - P*P)*(R - P*P);
		var complex_part = Math.sqrt(Math.abs(innerSquare));
		var result1 = Math.pow(Q*Q + complex_part*complex_part, 1/6)*Math.cos(1/3*(Math.atan2(complex_part, Q)+4*Math.PI));
		//var result2 = Math.pow(Q*Q + complex_part*complex_part, 1/6)*Math.cos(1/3*(Math.atan2(-complex_part, Q)+2*Math.PI));
		var result = 2*result1 + P;
		if (Math.abs(result) < 1e-10) {
			result = 0;
		}
		if (Math.abs(result-1) < 1e-10) {
			result = 1;
		}
		if (result > 1 || result < 0) {
			console.log("root choice error in constrained MLE!");
			console.log(result);
			console.log("S_x:"+S_x)
			console.log("S_y:"+S_y)
			console.log("n:"+n)
			console.log("d:"+d)
		}
		return result;
	}

	var bernoulli_twosided_alpha = function(b0, b1, indiff, simulateResult, samples) {
		var p1 = 0.5;
		var p2 = 0.5;
		if (!samples) samples = 10000;
		// calculate alpha error via importance sampling
		var alphas = []
		for (var i = 0;i < samples;i++) {
			var beta_alpha = 5;
			var beta_beta = 5;
			var p1_ran = jStat.jStat.beta.sample(beta_alpha,beta_beta);
			var p2_ran = jStat.jStat.beta.sample(beta_alpha,beta_beta);

			var res = simulateResult(p1_ran,p2_ran,b0,b1);
			if (res[0] == 'false') {
				var stoppingTime = res[3];
				var sum_x = res[1];
				var sum_y = res[2];
				var weight = Math.exp( logOp(sum_x, p2) + logOp(stoppingTime-sum_x, 1-p2) + jStat.jStat.betaln(beta_alpha, beta_beta) - jStat.jStat.betaln(beta_alpha+sum_x, beta_beta+stoppingTime-sum_x) + logOp(sum_y, p1) + logOp(stoppingTime-sum_y, 1-p1) + jStat.jStat.betaln(beta_alpha, beta_beta) - jStat.jStat.betaln(beta_alpha+sum_y, beta_beta+stoppingTime-sum_y) );
				alphas.push(weight);
			} else {
				alphas.push(0);
			}
		}
		return alphas;
	}

	var bernoulli_twosided_beta = function(b0, b1, indiff, simulateResult, samples) {
		if (!samples) samples = 10000;
		var betas = [];
		for (var i = 0;i < samples;i++) {
			var res = simulateResult(0,indiff,b0,b1);
			if (res[0] == 'true') {
				betas.push(1);
			} else {
				betas.push(0);
			}
		}
		return betas;
	}

	var bernoulli_twosided_LR_H0 = function(S_x, S_y, n, indiff) {
		var equal_mle = (S_x+S_y)/(2*n);
		// calculate unconstrained MLE, i.e. p1 and p2 can be unequal 
		var unc_mle_x = S_x/n;
		var unc_mle_y = S_y/n;

		var likRatio = Math.exp( (logOp(S_x,unc_mle_x) + logOp(n-S_x,1-unc_mle_x) + logOp(S_y,unc_mle_y) + logOp(n-S_y,1-unc_mle_y)) - (logOp(S_x,equal_mle) + logOp(n-S_x,1-equal_mle) + logOp(S_y,equal_mle) + logOp(n-S_y,1-equal_mle)));
		return likRatio;
	}

	var bernoulli_twosided_LR_HA = function(S_x, S_y, n, indiff) {
		var unc_mle_x = S_x/n;
		var unc_mle_y = S_y/n;

		if (Math.abs(unc_mle_x-unc_mle_y) > indiff) {
			return 1;
		}
		
		// find mle of p1 with constrain that |p1-p2| = d
		var pos = solveConstrainedBinomialMLE(S_x, S_y, n, indiff); // solves MLE of p1 with the constraint that p1 = p2 - d
		var neg = solveConstrainedBinomialMLE(S_x, S_y, n, -indiff); // solves MLE of p1 with the constraint that p1 = p2 + d

		var A_pos = roundToZero(pos);
		var B_pos = roundToZero(1-pos);
		var C_pos = roundToZero(pos + indiff);
		var D_pos = roundToZero(1-pos-indiff);
		var pos_llik = logOp(S_x,A_pos) + logOp(n-S_x,B_pos) + logOp(S_y,C_pos) + logOp(n-S_y,D_pos);

		var A_neg = roundToZero(neg);
		var B_neg = roundToZero(1-neg);
		var C_neg = roundToZero(neg - indiff);
		var D_neg = roundToZero(1-neg+indiff);
		var neg_llik = logOp(S_x,A_neg) + logOp(n-S_x,B_neg) + logOp(S_y,C_neg) + logOp(n-S_y,D_neg);
		
		if (pos_llik > neg_llik) {
			return Math.exp( logOp(S_x,unc_mle_x) + logOp(n-S_x,1-unc_mle_x) + logOp(S_y,unc_mle_y) + logOp(n-S_y,1-unc_mle_y) - pos_llik );
		} else {
			return Math.exp( logOp(S_x,unc_mle_x) + logOp(n-S_x,1-unc_mle_x) + logOp(S_y,unc_mle_y) + logOp(n-S_y,1-unc_mle_y) - neg_llik );
		}
	}

	var bernoulli_twosided_maxSamplesize = function(b0, b1, indiff) {
		var returnFunction = function() {
			// TODO : how to get threshold?
			var crossed = false;
			var L_na_thresholds = [];
			var L_nb_thresholds = [];
			var maxSample = 0;			
			for (var i = 0;!crossed;i++) {
				// start with S_y at Math.floor(0.5*n) and adjust S_y up until L_na crosses threshold (if it happens)
				var S_x = Math.floor(0.5*i);
				var S_y = Math.floor(0.5*i);
				var j = 0;
				while (S_y <= i && S_x >= 0) {
					if (bernoulli_twosided_LR_H0(S_x, S_y, i) >= b0) {
						L_na_thresholds[i] = Math.abs(S_x/i - S_y/i);
						break;
					}
					if (j % 2 == 0) S_y += 1;
					else S_x -= 1;
					j += 1;
				}
				// start with S_y at n and adjust S_Y down towards Math.floor(0.5*n) until L_nb crosses threshold (if it happens)
				var S_x = 0;
				var S_y = i;
				var j = 0;
				while (S_y >= Math.floor(0.5*i) && S_x <= Math.floor(0.5*i)) {
					if (bernoulli_twosided_LR_HA(S_x, S_y, i, indiff) >= b1) {
						L_nb_thresholds[i] = Math.abs(S_x/i - S_y/i);
						break;
					}
					if (j % 2 == 0) S_y -= 1;
					else S_x += 1;
					j += 1;
				}
				// if these crosses then we've reached worst case samplesize, so stop
				if (L_na_thresholds[i] <= L_nb_thresholds[i]) {
					maxSample = i;
					crossed = true;
				}
			}
			// write to file
			/*var fs = require('fs');
			var str1 = "c(";
			var str2 = "c(";
			for (var i = 0;i < maxSample;i++) {
				if (typeof L_na_thresholds[i] == 'undefined') {
					str1 += "NA,"
				} else {
					str1 += L_na_thresholds[i].toFixed(3)+","
				}
				if (typeof L_nb_thresholds[i] == 'undefined') {
					str2 += "NA,"
				} else {
					str2 += L_nb_thresholds[i].toFixed(3)+","
				}
			}
			fs.writeFile("./test.txt",str1+"),\n"+str2+")\n", function(err){});
			*/
			//return [maxSample, L_na_thresholds, L_nb_thresholds];

			return maxSample;
		}
		return returnFunction;
	}

	var bernoulli_twosided_simulateH0 = function(simRes, indiff, b0, b1) {
		var returnFun = function() {
			var res = simRes(0.5,0.5,b0,b1)[4];
			return res;
		}
		return returnFun;
	}

	var bernoulli_onesided_LR_H0 = function(S_x, S_y, n, indiff) {
		// nb! H0 is that p1 <= p2
		var unc_mle_x = S_x/n;
		var unc_mle_y = S_y/n;

		if (unc_mle_x-unc_mle_y <= -indiff) {
			return 1;
		}
		
		// p1 = p2 - indiff
		var pos = solveConstrainedBinomialMLE(S_x, S_y, n, indiff); // solves MLE of p1 with the constraint that p1 = p2 - d, i.e. p1 <= p2 - d

		var A_pos = roundToZero(pos);
		var B_pos = roundToZero(1-pos);
		var C_pos = roundToZero(pos + indiff);
		var D_pos = roundToZero(1-pos-indiff);
		var pos_llik = logOp(S_x,A_pos) + logOp(n-S_x,B_pos) + logOp(S_y,C_pos) + logOp(n-S_y,D_pos);

		return Math.exp( logOp(S_x,unc_mle_x) + logOp(n-S_x,1-unc_mle_x) + logOp(S_y,unc_mle_y) + logOp(n-S_y,1-unc_mle_y) - pos_llik );
	}

	var bernoulli_onesided_LR_HA = function(S_x, S_y, n, indiff) {
		// nb! HA is that p1 >= p2

		var unc_mle_x = S_x/n;
		var unc_mle_y = S_y/n;

		if (unc_mle_x-unc_mle_y >= indiff) {
			return 1;
		}
		
		// p1 = p2 + indiff
		var neg = solveConstrainedBinomialMLE(S_x, S_y, n, -indiff);

		var A_neg = roundToZero(neg);
		var B_neg = roundToZero(1-neg);
		var C_neg = roundToZero(neg - indiff);
		var D_neg = roundToZero(1-neg+indiff);
		var neg_llik = logOp(S_x,A_neg) + logOp(n-S_x,B_neg) + logOp(S_y,C_neg) + logOp(n-S_y,D_neg);

		return Math.exp( logOp(S_x,unc_mle_x) + logOp(n-S_x,1-unc_mle_x) + logOp(S_y,unc_mle_y) + logOp(n-S_y,1-unc_mle_y) - neg_llik );
	}

	var bernoulli_onesided_alpha = function(b0, b1, indiff, simulateResult, samples) {
		if (!samples) samples = 10000;
		var alphas = [];
		for (var i = 0;i < samples;i++) {
			var res = simulateResult(0.5-(indiff/2),0.5+(indiff/2),b0,b1);
			//var res = simulateResult(0,indiff/2,b0,b1);
			if (res[0] == 'false') {
				alphas.push(1);
			} else {
				alphas.push(0);
			}
		}
		return alphas;
	}

	var bernoulli_onesided_beta = function(b0, b1, indiff, simulateResult, samples) {
		if (!samples) samples = 10000;
		var betas = [];
		for (var i = 0;i < samples;i++) {
			var res = simulateResult(0.5+(indiff/2),0.5-(indiff/2),b0,b1);
			//var res = simulateResult(1,1-indiff,b0,b1);
			if (res[0] == 'true') {
				betas.push(1);
			} else {
				betas.push(0);
			}
		}
		return betas;
	}

	var bernoulli_onesided_maxSamplesize = function(b0, b1, indiff) {
		var returnFunction = function() {
			var crossed = false;
			var L_na_thresholds = [];
			var L_nb_thresholds = [];
			var maxSample = 0;			
			for (var i = 0;!crossed;i++) {
				// start with S_y at Math.floor(0.5*n) and adjust S_y up until L_na crosses threshold (if it happens)
				var S_x = Math.floor(0.5*i);
				var S_y = Math.floor(0.5*i);
				var j = 0;
				while (S_y >= 0 && S_x <= i) {
					if (onesided_LR_H0(S_x, S_y, i, indiff) >= b0) {
						L_na_thresholds[i] = S_x/i - S_y/i;
						break;
					}
					if (j % 2 == 0) S_y -= 1;
					else S_x += 1;
					j += 1;
				}
				// start with S_y at Math.floor(0.5*n) and adjust S_Y down until L_nb crosses threshold (if it happens)
				var S_x = Math.floor(0.5*i);
				var S_y = Math.floor(0.5*i);
				var j = 0;
				while (S_y <= i && S_x >= 0) {
					if (onesided_LR_HA(S_x, S_y, i, indiff) >= b1) {
						L_nb_thresholds[i] = S_x/i - S_y/i;
						break;
					}
					if (j % 2 == 0) S_y += 1;
					else S_x -= 1;
					j += 1;
				}
				// if these crosses then we've reached worst case samplesize, so stop
				if (L_na_thresholds[i] <= L_nb_thresholds[i]) {
					maxSample = i;
					crossed = true;
				}
			}
			// write to file
			var fs = require('fs');
			var str1 = "c(";
			var str2 = "c(";
			for (var i = 0;i < maxSample;i++) {
				if (typeof L_na_thresholds[i] == 'undefined') {
					str1 += "NA,"
				} else {
					str1 += L_na_thresholds[i].toFixed(3)+","
				}
				if (typeof L_nb_thresholds[i] == 'undefined') {
					str2 += "NA,"
				} else {
					str2 += L_nb_thresholds[i].toFixed(3)+","
				}
			}
			fs.writeFile("./test.txt",str1+"),\n"+str2+")\n", function(err){});

			return [maxSample, L_na_thresholds, L_nb_thresholds];
		}
		return returnFunction;
	}

	var bernoulli_onesided_simulateH0 = function(simRes, indiff, b0, b1) {
		var returnFun = function() {
			var res = simRes(0.5-indiff/2,0.5+indiff/2,b0,b1)[4];
			return res;
		}
		return returnFun;
	}

	functions['bernoulli'] = {
		'one-sided' : {
			'l_an' : bernoulli_onesided_LR_H0,
			'l_bn' : bernoulli_onesided_LR_HA,
			'alpha' : bernoulli_onesided_alpha,
			'beta' : bernoulli_onesided_beta,
			'max_samplesize' : bernoulli_onesided_maxSamplesize,
			'simulateH0' : bernoulli_onesided_simulateH0,
		},
		'two-sided' : {
			'l_an' : bernoulli_twosided_LR_H0,
			'l_bn' : bernoulli_twosided_LR_HA,
			'alpha' : bernoulli_twosided_alpha,
			'beta' : bernoulli_twosided_beta,
			'max_samplesize' : bernoulli_twosided_maxSamplesize,
			'simulateH0' : bernoulli_twosided_simulateH0,
		}
	}
	  
	thresholds['bernoulli'] = {
		'two-sided' : {
			0.05 : {
				0.05 : {
					0.1 : [140.5, 31],// approximate
				},
				0.10 : {
					0.4 : [68.6, 12.9], 
					0.2 : [98, 14.6], 
					0.1 : [139, 14.5],
					0.05 : [172, 15.5], 
					0.025 : [220, 15.5], 
					0.01 : [255, 15.7]
				},
				0.20 : {
					0.4 : [68, 5.4],
					0.2 : [90.5, 6.5],
					0.1 : [134, 6.9],
					0.05 : [168, 7.3],
					0.025 : [214, 7.4],
					0.01 : [254, 7.5]
				}
			}
		},
		'one-sided' : {
			0.05 : {
				0.05 : {
					0.2 : [39.1, 39.1],
					0.1 : [42, 42],
					0.05 : [70, 70],
					0.025 : [77, 77],
					0.01 : [95,95]
				}
			}
		}
	}

	tests['bernoulli'] = bernoulli_test;