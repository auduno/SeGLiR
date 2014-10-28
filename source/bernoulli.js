	
	/*** test for bernoulli proportions ***/

	var bernoulli_test = function(sides, indifference, type1_error, type2_error) {

		var b0, b1, stoppingTime;

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
		// this needs to be exchanged
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
			
			var result = checkTest(S_x, S_y, n, indiff, b0, b1);
			if (result) {
				finished = true;
				stoppingTime = n;
				L_an = result[1];
				return result[0];
			}
		}

		// get expected samplesize for some parameters
		/*this.expectedSamplesize = function(p1, p2, samples) {
			// simulate it enough times
			if (!samples) samples = 10000;
			console.log("calculating expected samplesize via simulation");
			var times = [];
			for (var i = 0;i < samples;i++) {
				var res = simulateResult(p1,p2,b0,b1)
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
				var res = simulateResult(p1,p2,b0,b1)
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
			
			// TODO : should I check for when both L_an and L_bn pass thresholds?

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
		this.LikH0 = LikH0;
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
		this.beta = beta;
		this.alpha = alpha;
		
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
		this.simulateResult = simulateResult;

		// initialization:
		  // calculate thresholds (unless they are stored in table)
		if (sides in thresholds['bernoulli'] && alpha_value in thresholds['bernoulli'][sides] && beta_value in thresholds['bernoulli'][sides][alpha_value] && indifference in thresholds['bernoulli'][sides][alpha_value][beta_value]) {
			b0 = thresholds['bernoulli'][sides][alpha_value][beta_value][indifference][0];
			b1 = thresholds['bernoulli'][sides][alpha_value][beta_value][indifference][1];
		} else {
			// calculate thresholds
			console.log("calculating thresholds via simulation")
			//var thr = optimize2d([alpha_value, beta_value], boundaryFun(indifference), [50,10], 0.001, 46000, 400000, 6, 1)
			//var thr = optimize2d([alpha_value, beta_value], boundaryFun(indifference), [98,14.5], 0.001, 46000, 1500000, 6, 1)
			var thr = optimize2d([alpha_value, beta_value], boundaryFun(indifference), [90,90], 0.001, 46000, 1500000, 6, 1)
			b0 = thr[0];
			b1 = thr[1];
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

		this.comparePower = function(p1,p2,samples) {
			var fixedPower = [];
			var seqPower = [];
			var samplesize = 52735;
			for (var i = 0;i < samples;i++) {
				// generate sample
				var res = simulateResult(p1,p2,b0,b1);
				if (res[0] == "true") {
					seqPower[i] = 1;
				} else {
					seqPower[i] = 0;
				}
				// test if reject H0
				var mle1 = 0;
				var mle2 = 0;
				for (var j = 0;j < samplesize;j++) {
					mle1 += generate(p1);
					mle2 += generate(p2);
				}
				mle1 /= samplesize;
				mle2 /= samplesize;
				var mlecomm = 0.5*mle1 + 0.5*mle2;
				var z_val = (mle1-mle2)/Math.sqrt(mlecomm*(1-mlecomm)*(1/samplesize + 1/samplesize))
				if (Math.abs(z_val) > 1.96) {
					fixedPower[i] = 0;
				} else {
					fixedPower[i] = 1;
				}
				// add
			}
			return [mean(seqPower), mean(fixedPower)];
		}

		this.compareErrors = function(p1, p2, samples) {
			var fixedErrors = 0;
			var seqErrors = 0;
			var allOver = 0;
			var samplesize = 52735;
			for (var i = 0;i < samples;i++) {
				var finished = false;
				var time = 0;
				var S_x = 0;
				var S_y = 0;
				var result;
				var S_x_fixed = 0;
				var S_y_fixed = 0;
				while (!finished) {
					S_x += generate(p1);
					S_y += generate(p2);
					time += 1;
					// test it
					var result = checkTest(S_x, S_y, time, indiff, b0, b1);
					if (result) finished = true;
					if (time == samplesize) {
						S_x_fixed = S_x;
						S_y_fixed = S_y;
					}
				}
				if (time > samplesize) {
					if (result[0] == "true") seqErrors += 1;
					allOver += 1;
					var mle1 = S_x_fixed/samplesize;
					var mle2 = S_y_fixed/samplesize;
					var mlecomm = 0.5*mle1 + 0.5*mle2;
					var z_val = (mle1-mle2)/Math.sqrt(mlecomm*(1-mlecomm)*(1/samplesize + 1/samplesize))
					if (Math.abs(z_val) <= 1.96) fixedErrors += 1;
				}
			}
			return [seqErrors/allOver, fixedErrors/allOver, allOver];
		}

		this.compareMAB = function(samples, length) {
			var glrRegrets = [];
			var thomRegrets = [];
			var classRegrets = [];
			for (var j = 0;j < length;j++) {
				glrRegrets[j] = 0;
				thomRegrets[j] = 0;
				classRegrets[j] = 0;
			}
			for (var i = 0;i < samples;i++) {
				// generate random p1 and p2
				var p1 = Math.random();
				var p2 = Math.random();
				
				// pull random center of alpha,beta
				/*var add =jStat.jStat.normal.sample(0,0.05);
				var p2 = p1 + add
				if (p2 > 1) {
					p2 = p1 - add
				}
				console.log(p2-p1)*/

				//console.log("p1:"+p1+",p2:"+p2);
				var prior1 = [1,1];
				var prior2 = [1,1];
				var S_x = 0;
				var S_y = 0;
				var time = 0;
				var finished = false;
				var choice = undefined;
				var correct_choice = p1 <= p2 ? 0 : 1;
				var glrRegret = 0;
				var thomRegret = 0;
				var classRegret = 0;
				for (var j = 0;j < length;j++) {
					var a = Math.random() < p1 ? 1 : 0;
					var b = Math.random() < p2 ? 1 : 0;
					// glr : pull alternate samples until finished (use one-sided with very very small indifference zone, alpha = ?, beta = ?)
					if (!finished) {
						if (j % 2 == 0) {
							S_x += a;
							/*if (time > 0) {
								var result = checkTest(S_x, S_y, time, indiff, b0, b1);
								if (result) {
									finished = true;
									if (result[0] == "true") choice = 0; // i.e. choose that p1 < p2
									else choice = 1; // i.e. choose that p1 > p2
								}
							}*/
							if (correct_choice == 0) {
								glrRegret += (p2-p1);
							}
							glrRegrets[j] += glrRegret;
						} else {
							S_y += b;
							time += 1;
							var result = checkTest(S_x, S_y, time, indiff, b0, b1);
							if (result) {
								finished = true;
								if (result[0] == "true") choice = 0; // i.e. choose that p1 < p2
								else choice = 1; // i.e. choose that p1 > p2
							}
							if (correct_choice == 1) {
								glrRegret += (p1-p2);
							}
							glrRegrets[j] += glrRegret;
						}
					} else {
						// do pull from choice
						if (correct_choice != choice) {
							if (correct_choice == 0) glrRegret += (p2-p1);
							if (correct_choice == 1) glrRegret += (p1-p2);
						}
						glrRegrets[j] += glrRegret;
					}
					// bayes bandit : pull from prior, choose, etc.
					var a_post_sample = jStat.jStat.beta.sample(prior1[0],prior1[1]);
					var b_post_sample = jStat.jStat.beta.sample(prior2[0],prior2[1]);
					if (a_post_sample >= b_post_sample) {
						// choose sample a
						prior1[0] += a;
						prior1[1] += (1-a);
						if (correct_choice == 0) {
							thomRegret += (p2-p1);
						}
						thomRegrets[j] += thomRegret;
					} else {
						// choose sample b
						prior2[0] += b;
						prior2[1] += (1-b);
						if (correct_choice == 1) {
							thomRegret += (p1-p2);
						}
						thomRegrets[j] += thomRegret;
					}
					if (j % 2 == 0) {
						if (correct_choice == 0) {
							classRegret += (p2-p1);
						} else {
							classRegret += (p1-p2);
						}
					}
					classRegrets[j] += classRegret;
				}
			}
			for (var i = 0;i < length;i++) {
				glrRegrets[i] /= samples;
				thomRegrets[i] /= samples;
				classRegrets[i] /= samples;
			}
			// return array of mean regret
			return [glrRegrets, thomRegrets, classRegrets];
		}

		this.compareMAB2 = function(samples, length,eps) {
			var glrRegrets = [];
			var thomRegrets = [];
			for (var j = 0;j < length;j++) {
				glrRegrets[j] = 0;
				thomRegrets[j] = 0;
			}
			for (var i = 0;i < samples;i++) {
				// generate random p1 and p2
				var p1 = Math.random();
				var p2 = Math.random();
				//console.log("p1:"+p1+",p2:"+p2);
				var prior1 = [1,1];
				var prior2 = [1,1];
				var S_x = 0;
				var S_y = 0;
				var time = 0;
				var finished = false;
				var choice = undefined;
				var correct_choice = p1 <= p2 ? 0 : 1;
				var glrRegret = 0;
				var thomRegret = 0;
				for (var j = 0;j < length;j++) {
					var a = Math.random() < p1 ? 1 : 0;
					var b = Math.random() < p2 ? 1 : 0;
					// glr : pull alternate samples until finished (use one-sided with very very small indifference zone, alpha = ?, beta = ?)
					if (!finished) {
						if (j % 2 == 0) {
							S_x += a;
							if (time > 0) {
								var L_an = LikH0(S_x, S_y, time, indiff);
								if (L_an > Math.log(time*2 + 1)/eps) {
									finished = true;
									if (S_x/(time+1) < S_y/time) choice = 0; // i.e. choose that p1 < p2
									else choice = 1; // i.e. choose that p1 > p2
								}
							}
							if (correct_choice == 0) {
								glrRegret += (p2-p1);
							}
							glrRegrets[j] += glrRegret;
						} else {
							S_y += b;
							time += 1;
							var L_an = LikH0(S_x, S_y, time, indiff);
							if (L_an > Math.log(time*2)/eps) {
								finished = true;
								if (S_x/time < S_y/time) choice = 0; // i.e. choose that p1 < p2
								else choice = 1; // i.e. choose that p1 > p2
							}
							if (correct_choice == 1) {
								glrRegret += (p1-p2);
							}
							glrRegrets[j] += glrRegret;
						}
					} else {
						// do pull from choice
						if (correct_choice != choice) {
							if (correct_choice == 0) glrRegret += (p2-p1);
							if (correct_choice == 1) glrRegret += (p1-p2);
						}
						glrRegrets[j] += glrRegret;
					}
					// bayes bandit : pull from prior, choose, etc.
					var a_post_sample = jStat.jStat.beta.sample(prior1[0],prior1[1]);
					var b_post_sample = jStat.jStat.beta.sample(prior2[0],prior2[1]);
					if (a_post_sample >= b_post_sample) {
						// choose sample a
						prior1[0] += a;
						prior1[1] += (1-a);
						if (correct_choice == 0) {
							thomRegret += (p2-p1);
						}
						thomRegrets[j] += thomRegret;
					} else {
						// choose sample b
						prior2[0] += b;
						prior2[1] += (1-b);
						if (correct_choice == 1) {
							thomRegret += (p1-p2);
						}
						thomRegrets[j] += thomRegret;
					}
				}
			}
			for (var i = 0;i < length;i++) {
				glrRegrets[i] /= samples;
				thomRegrets[i] /= samples;
			}
			// return array of mean regret
			return [glrRegrets, thomRegrets];
		}

		this.emilie = function(p1, p2, eps, samples) {
			// simulate it enough times
			if (!samples) samples = 10000;
			var times = [];
			for (var i = 0;i < samples;i++) {
				var S_x = 0;
				var S_y = 0;
				var time = 0;
				var finished = false;
				while (!finished) {
					S_x += generate(p1);
					S_y += generate(p2);
					time += 1;
					var Lna = LikH0(S_x, S_y, time);
					if (Lna > Math.log(2*time)/eps) {
						finished = true;
					}
				}
				times.push(time);
			}
			times.sort(function(a,b){return a-b});
			return [times[samples*0.05], mean(times), times[samples*0.95]];
		}
	}

	// private functions

	var solveConstrainedBinomialMLE = function(S_x, S_y, n, d) {
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
		// TODO : should we include std.dev.?
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
		// TODO : should we include std.dev.?
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
		
		var pos = solveConstrainedBinomialMLE(S_x, S_y, n, indiff);
		var neg = solveConstrainedBinomialMLE(S_x, S_y, n, -indiff);

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

		if (pos_llik < neg_llik) {
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
					if (twosided_LR_H0(S_x, S_y, i) >= b0) {
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
					if (twosided_LR_HA(S_x, S_y, i, indiff) >= b1) {
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

			return [maxSample, L_na_thresholds, L_nb_thresholds];
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
		var pos = solveConstrainedBinomialMLE(S_x, S_y, n, indiff);

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
		// TODO : should we include std.dev.?
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
		// TODO : should we include std.dev.?
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
				0.10 : {
					0.4 : [68.6, 12.9], // check
					0.2 : [98, 14.6], // check
					0.1 : [139, 14.5], // check
					0.05 : [172, 15.5], // check
					0.025 : [220, 15.5], // check
					0.01 : [255, 15.7], // check
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