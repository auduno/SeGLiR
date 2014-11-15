var jStat = require('jStat');

// main class
var glr = function() {
	var functions = {}

	// precalculated thresholds
	var thresholds = {}

	var tests = {}

	this.test = function(type) {
		if (type in tests) {
			return new tests[type](arguments[1], arguments[2], arguments[3], arguments[4], arguments[5], arguments[6], arguments[7]);
		} else {
			console.log("No test of type '"+type+"'.");
		}
	}

	/** generic utility functions **/

	/*
	 * use gradient descent in 2d to find parameters x,y so that fun(x,y) == vec
	 */
	var optimize2d = function(vec, fun, init_points, epsilon, nsamples, max_samples, gradient_evaluation_width, lower_limit, upper_limit, verbose) {
		if (typeof(verbose) === 'undefined') {
			verbose = true;
		}

		// clone variables
		vec = [vec[0],vec[1]];
		init_points = [init_points[0],init_points[1]];

		var gradient, point, next_point;
		var samples = 0;
		var diff = Infinity;
		var est_point = init_points;
		if (verbose) console.log(est_point);
		while (samples < max_samples && diff > epsilon) {
			// evaluate at current point

			// TODO : need to fix if est_points are close to boundaries, closer than gradient evaluation width

			// estimate gradient
			var l_point = [est_point[0] - gradient_evaluation_width/2, est_point[1]];
			var r_point = [est_point[0] + gradient_evaluation_width/2, est_point[1]];
			var d_point = [est_point[0], est_point[1] - gradient_evaluation_width/2];
			var u_point = [est_point[0], est_point[1] + gradient_evaluation_width/2];

			var enoughSamples = false;
			var l_samples = [[],[]];
			var r_samples = [[],[]];
			var u_samples = [[],[]];
			var d_samples = [[],[]];

			// check whether p-values are within epsilon from true p-value
			// i.e. sample until we can say with some certainty whether they are or not
			// if they are, stop estimation
			var curval = [[],[]];
			var withinZero = false;
			var curpoints;
			for (var i = 0;!withinZero && samples < max_samples;i++) {
				num_samples = nsamples*Math.pow(2,i);
				var new_curval = fun(est_point, num_samples);
				curval[0] = curval[0].concat(new_curval[0]);
				curval[1] = curval[1].concat(new_curval[1]);
				
				curpoints = [mean(curval[0]), mean(curval[1])];
				if (verbose) console.log("alpha : "+curpoints[0]+" +- "+(4*std(curval[0])/Math.sqrt(curval[0].length)));
				if (verbose) console.log("beta : "+curpoints[1]+" +- "+(4*std(curval[1])/Math.sqrt(curval[1].length)));

				// if CI does not contain 0 OR CI is smaller than 2*epsilon, stop
				var ci_halfwidth_0 = (4*std(curval[0])/Math.sqrt(curval[0].length));
				var ci_halfwidth_1 = (4*std(curval[1])/Math.sqrt(curval[1].length));
				var lower0 = curpoints[0] - ci_halfwidth_0;
				var upper0 = curpoints[0] + ci_halfwidth_0;
				var lower1 = curpoints[1] - ci_halfwidth_1;
				var upper1 = curpoints[1] + ci_halfwidth_1;
				if (sign(lower0-vec[0]) === sign(upper0-vec[0]) || sign(lower1-vec[1]) === sign(upper1-vec[1])) {
					withinZero = true;
				} else if ( ci_halfwidth_0 < epsilon && ci_halfwidth_1 < epsilon ) {
					withinZero = true;
					diff = Math.sqrt( Math.pow(ci_halfwidth_0,2) + Math.pow(ci_halfwidth_1,2) );
				}
				samples += num_samples;

				if (verbose) console.log("checked current estimate, samples:"+samples)
			}
			if (diff < epsilon || samples > max_samples) {
				break;
			}

			var i = 0;
			while (!enoughSamples && samples < max_samples) {
				num_samples = nsamples*Math.pow(2,i);
				// get samples from points
				var new_l_samples = fun(l_point, num_samples);
				var new_r_samples = fun(r_point, num_samples);
				var new_u_samples = fun(u_point, num_samples);
				var new_d_samples = fun(d_point, num_samples);


				l_samples[0] = l_samples[0].concat(new_l_samples[0]);
				l_samples[1] = l_samples[1].concat(new_l_samples[1]);
				r_samples[0] = r_samples[0].concat(new_r_samples[0]);
				r_samples[1] = r_samples[1].concat(new_r_samples[1]);
				u_samples[0] = u_samples[0].concat(new_u_samples[0]);
				u_samples[1] = u_samples[1].concat(new_u_samples[1]);
				d_samples[0] = d_samples[0].concat(new_d_samples[0]);
				d_samples[1] = d_samples[1].concat(new_d_samples[1]);
				
				if (verbose) console.log("length samples : "+l_samples[0].length);

				samples += num_samples;
				if (verbose) console.log(samples);

				var l_0_mean = mean(l_samples[0]);
				var l_1_mean = mean(l_samples[1]);
				var r_0_mean = mean(r_samples[0]);
				var r_1_mean = mean(r_samples[1]);
				var u_0_mean = mean(u_samples[0]);
				var u_1_mean = mean(u_samples[1]);
				var d_0_mean = mean(d_samples[0]);
				var d_1_mean = mean(d_samples[1]);
				var b1p1_gradient_mean = (r_0_mean-l_0_mean)/gradient_evaluation_width;
				var b1p2_gradient_mean = (u_0_mean-d_0_mean)/gradient_evaluation_width;
				var b2p1_gradient_mean = (r_1_mean-l_1_mean)/gradient_evaluation_width;
				var b2p2_gradient_mean = (u_1_mean-d_1_mean)/gradient_evaluation_width;
				//console.log("gradient : "+gradient_mean+" +- "+(4*( (std(l_samples)+std(r_samples))/Math.pow(gradient_evaluation_width,2) )/Math.sqrt(nsamples)));
				
				/*console.log("b1p1 gradient : "+b1p1_gradient_mean+" +- "+(4*( (std(r_samples[0])+std(l_samples[0]))/Math.pow(gradient_evaluation_width,2) )/Math.sqrt(nsamples)));
				console.log("b1p2 gradient : "+b1p2_gradient_mean+" +- "+(4*( (std(u_samples[0])+std(d_samples[0]))/Math.pow(gradient_evaluation_width,2) )/Math.sqrt(nsamples)));
				console.log("b2p1 gradient : "+b2p1_gradient_mean+" +- "+(4*( (std(r_samples[1])+std(l_samples[1]))/Math.pow(gradient_evaluation_width,2) )/Math.sqrt(nsamples)));
				console.log("b2p2 gradient : "+b2p2_gradient_mean+" +- "+(4*( (std(u_samples[1])+std(d_samples[1]))/Math.pow(gradient_evaluation_width,2) )/Math.sqrt(nsamples)));*/

				/*var b1p1_cov = 0;
				for (var i = 0;i < r_samples[0].length;i++) {
					b1p1_cov += (r_samples[0][i]-r_0_mean)*(l_samples[0][i]-l_0_mean);
				}
				b1p1_cov /= (r_samples[0].length-1)
				console.log("b1p1_cov:"+b1p1_cov);
				console.log(Math.pow(std(r_samples[0]),2))
				console.log(Math.pow(std(l_samples[0]),2))
				console.log((gradient_evaluation_width*gradient_evaluation_width*r_samples[0].length))
				console.log(( Math.pow(std(r_samples[0]),2)+Math.pow(std(l_samples[0]),2) )/(gradient_evaluation_width*gradient_evaluation_width*r_samples[0].length))
				var b1p1_sd = ( Math.pow(std(r_samples[0]),2)+Math.pow(std(l_samples[0]),2) )/(gradient_evaluation_width*gradient_evaluation_width*r_samples[0].length) - (2*b1p1_cov)/(gradient_evaluation_width*gradient_evaluation_width);
				console.log((2*b1p1_cov)/(gradient_evaluation_width*gradient_evaluation_width))*/

				var b1p1_var = ( Math.pow(std(r_samples[0]),2)+Math.pow(std(l_samples[0]),2) )/(gradient_evaluation_width*gradient_evaluation_width*r_samples[0].length);
				var b1p2_var = ( Math.pow(std(u_samples[0]),2)+Math.pow(std(d_samples[0]),2) )/(gradient_evaluation_width*gradient_evaluation_width*r_samples[0].length);
				var b2p1_var = ( Math.pow(std(r_samples[1]),2)+Math.pow(std(l_samples[1]),2) )/(gradient_evaluation_width*gradient_evaluation_width*r_samples[0].length);
				var b2p2_var = ( Math.pow(std(u_samples[1]),2)+Math.pow(std(d_samples[1]),2) )/(gradient_evaluation_width*gradient_evaluation_width*r_samples[0].length);

				if (verbose) console.log("b1p1 gradient : "+b1p1_gradient_mean+" +- "+(4*Math.sqrt(b1p1_var)) );
				if (verbose) console.log("b1p2 gradient : "+b1p2_gradient_mean+" +- "+(4*Math.sqrt(b1p2_var)) );
				if (verbose) console.log("b2p1 gradient : "+b2p1_gradient_mean+" +- "+(4*Math.sqrt(b2p1_var)) );
				if (verbose) console.log("b2p2 gradient : "+b2p2_gradient_mean+" +- "+(4*Math.sqrt(b2p2_var)) );
				
				enoughSamples = true;
				if (sign(b1p1_gradient_mean+4*Math.sqrt(b1p1_var)) != sign(b1p1_gradient_mean-4*Math.sqrt(b1p1_var))) enoughSamples = false;
				if (sign(b2p2_gradient_mean+4*Math.sqrt(b2p2_var)) != sign(b2p2_gradient_mean-4*Math.sqrt(b2p2_var))) enoughSamples = false;
				//if (sign(b1p2_gradient_mean+4*Math.sqrt(b1p2_var)) != sign(b1p2_gradient_mean-4*Math.sqrt(b1p2_var))) enoughSamples = false;
				//if (sign(b2p1_gradient_mean+4*Math.sqrt(b2p1_var)) != sign(b2p1_gradient_mean-4*Math.sqrt(b2p1_var))) enoughSamples = false;
				
				if (verbose) console.log("b1p1*b2p2:"+(b1p1_gradient_mean*b2p2_gradient_mean));
				if (verbose) console.log("b1p2*b2p1:"+(b1p2_gradient_mean*b2p1_gradient_mean));
				if (verbose) console.log("b2p2*(v0-c0):"+( b2p2_gradient_mean*(vec[0]-curpoints[0]) ));
				if (verbose) console.log("b1p2*(v1-c1):"+( b1p2_gradient_mean*(vec[1]-curpoints[1]) ));

				i += 1;
				if (verbose) console.log("getting gradients, samples:"+samples)
			}

			// extrapolate where point lies with simple linear function
			var mult = 1/(b1p1_gradient_mean*b2p2_gradient_mean - b1p2_gradient_mean*b2p1_gradient_mean);
			next_point = [mult,mult];
			next_point[0] *= ( b2p2_gradient_mean*(vec[0]-curpoints[0]) - b1p2_gradient_mean*(vec[1]-curpoints[1]) );
			next_point[1] *= ( b1p1_gradient_mean*(vec[1]-curpoints[1]) - b2p1_gradient_mean*(vec[0]-curpoints[0]) );

			// calculate difference between new point and estimated point
			//diff = Math.sqrt(next_point[0]*next_point[0] + next_point[1]*next_point[1]);

			//next_point = est_point + 0.5*(next_point-est_point);
			if (verbose) console.log(next_point);
			est_point[0] += next_point[0];
			est_point[1] += next_point[1];

			if (upper_limit) {
				if (est_point[0] > upper_limit) {
					est_point[0] = upper_limit;
				}
				if (est_point[1] > upper_limit) {
					est_point[1] = upper_limit;
				}
			}
			if (lower_limit) {
				if (est_point[0] < lower_limit) {
					est_point[0] = lower_limit;
				}
				if (est_point[1] < lower_limit) {
					est_point[1] = lower_limit;
				}
			}

			if (verbose) console.log(est_point);
		}
		// calculate estimate of final value
		var curval = fun(est_point, nsamples);
		curpoints = [mean(curval[0]), mean(curval[1])];

		if (verbose) console.log("alpha : "+curpoints[0]+" +- "+(4*std(curval[0])/Math.sqrt(nsamples)));
		if (verbose) console.log("beta : "+curpoints[1]+" +- "+(4*std(curval[1])/Math.sqrt(nsamples)));

		return est_point;
	}

	var mean = function(seq) {
		var sum = 0;
		for (var i = 0;i < seq.length;i++) {
			sum += seq[i];
		}
		return sum/seq.length;
	}

	var boot_std = function(seq,n) {
		var sl = seq.length;
		var boot_means = [];
		for (var i = 0;i < n;i++) {
			var boot_seq = [];
			for (var j = 0;j < sl;j++) {
				var ind = Math.floor(Math.random()*sl);
				boot_seq.push(seq[ind]);
			}
			boot_means.push(mean(boot_seq));
		}
		return std(boot_means);
	}

	var std = function(seq) {
		var mean_seq = mean(seq);
		var sum = 0;
		for (var i = 0;i < seq.length;i++) {
			sum += (seq[i]-mean_seq)*(seq[i]-mean_seq);
		}
		sum /= seq.length;
		return Math.sqrt(sum);
	}

	var sign = function(x) {
	    if( +x === x ) {
	        return (x === 0) ? x : (x > 0) ? 1 : -1;
	    }
	    return NaN;
	}

	var roundToZero = function(x) {
		if (Math.abs(x) < 1e-10) {
			return 0;
		}
		return x;
	}

	var logOp = function(mult, logvar) {
		// by default 0*Infinite = NaN in javascript, so we make a custom operator where this will be equal to 0
		if (mult == 0) {
			return 0;
		}
		var logged = Math.log(logvar);
		if (!isFinite(logged)) {
			return logged;
		} 
		return mult*logged;
	}

	
	/*** test for bernoulli proportions ***/

	var bernoulli_test = function(sides, indifference, type1_error, type2_error) {

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

		var betaProb = function(prior1,prior2) {
			// from http://www.evanmiller.org/bayesian-ab-testing.html
			var a1 = prior1[0];
			var b1 = prior1[1];
			var a2 = prior2[0];
			var b2 = prior2[1];
			var sum = 0
			for (var i = 0;i < a2;i++) {
				//sum += jStat.jStat.betafn(i + a1, b1 + b2)/( (b2+i)*jStat.jStat.betafn(1 + i, b2) * jStat.jStat.betafn(a1, b1) );
				sum += Math.exp( jStat.jStat.betaln(i + a1, b1 + b2) - Math.log(b2+i) - jStat.jStat.betaln(1 + i, b2) - jStat.jStat.betaln(a1, b1) );
			}
			return 1-sum;
		}
		this.betaProb = betaProb;

		this.expectedThomError = function(p1, p2, samples) {
			var errors = 0;
			var correct_choice = p1 <= p2 ? 1 : 0;
			for (var i = 0;i < samples;i++) {
				//console.log("p1:"+p1+",p2:"+p2);
				var prior1 = [1,1];
				var prior2 = [1,1];
				var S_x = 0;
				var S_y = 0;
				var time = 0;
				var finished = false;
				var choice = undefined;
				var j = 0;
				while (!finished) {
					j += 1;
					var a = Math.random() < p1 ? 1 : 0;
					var b = Math.random() < p2 ? 1 : 0;
					// bayes bandit : pull from prior, choose, etc.
					var a_post_sample = jStat.jStat.beta.sample(prior1[0], prior1[1]);
					var b_post_sample = jStat.jStat.beta.sample(prior2[0], prior2[1]);
					if (a_post_sample >= b_post_sample) {
						// choose sample a
						prior1[0] += a;
						prior1[1] += (1-a);
					} else {
						// choose sample b
						prior2[0] += b;
						prior2[1] += (1-b);
					}
					var prob = betaProb(prior1,prior2)
					if (prob > 0.95) {
						finished = true;
						choice = 0;
					} else if (1-prob > 0.95) {
						finished = true;
						choice = 1;
					}
				}
				if (choice != correct_choice) {
					errors += 1;
				}
				console.log(i+" : "+(errors/(i+1)));
			}
			return errors/samples;
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
					0.4 : [68.6, 12.9], 
					0.2 : [98, 14.6], 
					0.1 : [139, 14.5],
					0.05 : [172, 15.5], 
					0.025 : [220, 15.5], 
					0.01 : [255, 15.7]
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
			var lower_est_bc = optimize2d(lower_est, biasFun(), lower_est, 0.005, 16400, 590000, 0.02, undefined, undefined, false);
			var upper_est_bc = optimize2d(upper_est, biasFun(), upper_est, 0.005, 16400, 590000, 0.02, undefined, undefined, false);

			return [(lower_est_bc[0]-lower_est_bc[1]),(upper_est_bc[0]-upper_est_bc[1])];
		}

		// get estimate (only when test is done)
		  // use bias-reduction
		this.estimate = function() {
			if (!finished) {
				return undefined;
			}
			var ests = optimize2d([S_x/n, S_y/n], biasFun(), [S_x/n, S_y/n], 0.005, 16400, 590000, 0.02, undefined, undefined, false);
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
				if (typeof points['x'] === 'number') x_data.push(points['x']);
				if (typeof points['y'] === 'number') y_data.push(points['y']);
			} else {
				if (typeof points['x'] === 'number' && typeof points['y'] === 'number') {
					if (x_data.length == y_data.length) {
						S_x += points['x'];
						S_y += points['y'];
						S_x2 += points['x']*points['x'];
						S_y2 += points['y']*points['y'];
						n += 1;
					} else if (x_data.length > y_data.length) {
						S_x += x_data[n];
						S_y += points['y'];
						S_x2 += x_data[n]*x_data[n];
						S_y2 += points['y']*points['y'];
						n += 1;
					} else {
						S_x += points['x'];
						S_y += y_data[n];
						S_x2 += points['x']*points['x'];
						S_y2 += y_data[n]*y_data[n];
						n += 1;
					}
					x_data.push(points['x'])
					y_data.push(points['y'])
				} else if (typeof points['x'] === 'number') {
					if (x_data.length < y_data.length) {
						S_x += points['x'];
						S_y += y_data[n];
						S_x2 += points['x']*points['x'];
						S_y2 += y_data[n]*y_data[n];
						n += 1;
					}
					x_data.push(points['x']);
				} else if (typeof points['y'] === 'number') {
					if (x_data.length > y_data.length) {
						S_x += x_data[n];
						S_y += points['y'];
						S_x2 += x_data[n]*x_data[n];
						S_y2 += points['y']*points['y'];
						n += 1;
					}
					y_data.push(points['y']);
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
			var est_var;
			if (typeof(variance) == "undefined") {
				est_var = (n*S_x2 - S_x*S_x + n*S_y2 - S_y*S_y)/(2*n*(n-1));
			} else {
				est_var = variance;
			}

			var outfun = function(parameters, n) {
				var results_p1 = []
				var results_p2 = []
				for (var i = 0;i < n;i++) {
					// generate sequences
					var res = simulateResult([parameters[0],est_var], [parameters[1],est_var], b0, b1);
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
			var thr = optimize2d([alpha_value, beta_value], boundaryFun(indifference), [100,100], 0.001, 46000, 1500000, 6, 1)
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
		var starttime = (new Date()).getTime();
		for (var i = 0;i < samples;i++) {
			var res = simulateResult([0,var_val],[0,var_val],b0,b1);
			if (res[0] == 'false') {
				alphas.push(1);
			} else {
				alphas.push(0);
			}
		}
		console.log("time:"+( (new Date()).getTime()-starttime ))
		console.log("mean:"+mean(alphas));
		console.log("std_err:"+boot_std(alphas,1000));
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

	var normal_uv_twosided_alpha_imp = function(b0, b1, indiff, var_val, simulateResult, samples) {
		if (!samples) samples = 10000;
		var alphas = [];
		var beta = 1; // precision/inverse-variance of the importance sampling distribution
		//var starttime = (new Date()).getTime();
		for (var i = 0;i < samples;i++) {
			var z = jStat.jStat.normal.sample(0,Math.sqrt(1/beta));

			var finished = false;
			var S_x = 0;
			var S_x2 = 0
			var n = 0;
			var result = undefined;
			while (!finished) {
				n += 1;
				if (n > 1) {
					// pull xs from N(0,2*var_val)
					var sample = jStat.jStat.normal.sample(z,Math.sqrt(2*var_val));
					S_x += sample;
					S_x2 += sample*sample;
					// test on simplified boundaries
					var L_na = Math.exp( n/2 * Math.log((n*S_x2)/(n*S_x2 - S_x*S_x)) );
					if (L_na >= b0) {
						finished = true;
						result = "false"
					}
					if (Math.abs(S_x/n) < indiff) {
						if (S_x/n > 0) {
							var L_nb = Math.exp( n/2 * Math.log( (n*(S_x2 - 2*indiff*S_x + n*indiff*indiff)) / (n*S_x2 - S_x*S_x) ))
						} else {
							var L_nb = Math.exp( n/2 * Math.log( (n*(S_x2 + 2*indiff*S_x + n*indiff*indiff)) / (n*S_x2 - S_x*S_x) ))
						}
						if (L_nb >= b1) {
							finished = true;
							result = "true";
						}
					}
				}
				/*if (n % 10 == 0 && n > 0) {
					console.log("****");
					console.log("S_x:"+S_x);
					console.log("S_x2:"+S_x2);
					console.log("n:"+n);
					console.log("L_na:"+L_na);
					console.log("L_nb:"+L_nb);
					break;
				}*/
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

	var normal_uv_twosided_alpha_imp2 = function(b0, b1, indiff, var_val, simulateResult, samples) {
		if (!samples) samples = 10000;
		var alphas = [];
		var beta = 10;
		var starttime = (new Date()).getTime();
		for (var i = 0;i < samples;i++) {
			var mu_1 = jStat.jStat.normal.sample(0,Math.sqrt(1/beta));
			var mu_2 = jStat.jStat.normal.sample(0,Math.sqrt(1/beta));
			var res = simulateResult([mu_1,var_val],[mu_2,var_val],b0,b1);
			if (res[0] == 'false') {
				var S_x = res[1];
				var S_y = res[2];
				var time = res[5];
				var bv = beta*var_val;
				
				//var weight = (bv/(time+bv))*Math.exp((1/var_val)*((S_x*S_x - 2*S_x*S_y + S_y*S_y)/(2*(time+bv))));
				var weight = (bv/(time+bv))*Math.exp((1/var_val)*((S_x*S_x + S_y*S_y)/(2*(time+bv))));
				//var weight = (bv/(time+bv))*Math.exp((1/var_val) * ( -(S_x + S_y)*(S_x + S_y)/(4*time) + (S_x*S_x + S_y*S_y)/(2*(time+bv)) ));
				
				//var weight = beta/(2*time + beta)*Math.exp( (S_x*S_x + S_y*S_y)/(time + 0.5*beta) - (S_x + S_y)*(S_x + S_y)/(2*time) );
				//var weight = Math.sqrt(beta/(time+beta))*Math.exp(time*time* ((S_x-S_y)/(2*time))*((S_x-S_y)/(2*time)) / (2*(time+beta)));
				alphas.push(1/weight);
			} else {
				alphas.push(0);
			}
		}
		console.log("time:"+( (new Date()).getTime()-starttime ))
		console.log("mean:"+mean(alphas));
		console.log("std_err:"+boot_std(alphas,1000));
		return alphas;
		// TODO : should we include std.dev.?
	}

	var normal_uv_twosided_alpha_imp3 = function(b0, b1, indiff, var_val, simulateResult, samples) {
		if (!samples) samples = 10000;
		var alphas = [];
		var beta = 1; // precision/inverse-variance of the importance sampling distribution
		var starttime = (new Date()).getTime();
		for (var i = 0;i < samples;i++) {
			var z = jStat.jStat.normal.sample(0,Math.sqrt(1/beta));
			var mu_1 = -z/2;
			var mu_2 = z/2;
			var res = simulateResult([mu_1,var_val],[mu_2,var_val],b0,b1);
			if (res[0] == 'false') {
				var S_x = res[1];
				var S_y = res[2];
				var time = res[5];
				var b2v = 2*beta*var_val;
				
				var weight = (Math.sqrt(b2v)/Math.sqrt(time + b2v))*Math.exp( (S_y-S_x)*(S_y-S_x)/(4*var_val*(time+b2v)) );
				alphas.push(1/weight);
			} else {
				alphas.push(0);
			}
		}
		console.log("time:"+( (new Date()).getTime()-starttime ))
		console.log("mean:"+mean(alphas));
		console.log("std_err:"+boot_std(alphas,1000));
		return alphas;
		// TODO : should we include std.dev.?
	}

	var normal_twosided_beta = function(b0, b1, indiff, var_val, simulateResult, samples) {
		if (!samples) samples = 10000;
		var betas = [];
		var starttime = (new Date()).getTime();
		for (var i = 0;i < samples;i++) {
			var res = simulateResult([-indiff/2,var_val],[indiff/2,var_val],b0,b1);
			if (res[0] == 'true') {
				betas.push(1);
			} else {
				betas.push(0);
			}
		}
		console.log("time:"+( (new Date()).getTime()-starttime ))
		console.log("mean:"+mean(betas));
		console.log("std_err:"+boot_std(betas,1000));
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

	var normal_uv_twosided_beta_imp2 = function(b0, b1, indiff, var_val, simulateResult, samples) {
		if (!samples) samples = 10000;
		var betas = [];
		var starttime = (new Date()).getTime();
		for (var i = 0;i < samples;i++) {	
			var res = simulateResult([0,var_val],[0,var_val],b0,b1);
			if (res[0] == 'true') {
				var S_x = res[1];
				var S_y = res[2];
				var time = res[5];
				//var weight = Math.exp((1/(var_val*2))*(indiff*S_y + indiff*S_x - time*indiff*indiff/2));
				var weight = Math.exp((1/(var_val*2))*(indiff*S_y - indiff*S_x - time*indiff*indiff/2));
				betas.push(weight);
			} else {
				betas.push(0);
			}
		}
		console.log("time:"+( (new Date()).getTime()-starttime ))
		console.log("mean:"+mean(betas));
		console.log("std_err:"+boot_std(betas,1000));
		return betas;
		// TODO : should we include std.dev.?
	}

	var normal_uv_twosided_LR_H0 = function(S_x, S_y, S_x2, S_y2, n, indiff) {
		if (n == 1) {
			// TODO : when n = 1, mle_llik is always -Infinity, so we have to enforce n > 1
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
			// TODO : when n = 1, mle_llik is always -Infinity, so we have to enforce n > 1
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

	var normal_onesided_alpha_imp = function(b0, b1, indiff, var_val, simulateResult, samples) {
		if (!samples) samples = 10000;
		var alphas = [];
		var beta = 1; // precision/inverse-variance of the importance sampling distribution
		//var starttime = (new Date()).getTime();
		for (var i = 0;i < samples;i++) {
			var z = jStat.jStat.normal.sample(0,Math.sqrt(1/beta));
			z = Math.abs(z);

			var finished = false;
			var S_x = 0;
			var n = 0;
			var result = undefined;
			while (!finished) {
				n += 1;
				// pull xs from N(0,2*var_val)
				S_x += jStat.jStat.normal.sample(z,Math.sqrt(2*var_val));

				var mle = S_x/n;
				// test on simplified boundaries
				if (mle < -indiff) {
					var L_na = 1;
				} else {
					var L_na = Math.exp( n*(S_x/n + indiff)*(S_x/n + indiff)/(4*var_val) );
				}
				if (L_na >= b0) {
					finished = true;
					result = "false"
				}
				if (mle > indiff) {
					var L_nb = 1;
				} else {
					var L_nb = Math.exp( n*(S_x/n - indiff)*(S_x/n - indiff)/(4*var_val) );
				}
				if (L_nb >= b1) {
					finished = true;
					result = "true";
				}
			}

			if (result == 'false') {
				var delta = (n + 2*var_val*beta)/(4*var_val);
				var pt1 = Math.sqrt(beta/(2*delta));
				var pt2 = Math.exp((2*indiff*S_x + n*indiff*indiff)/(4*var_val) + (S_x*S_x)/(16*delta*var_val*var_val));
				var pt3 = jStat.jStat.erf(S_x/(4*var_val*Math.sqrt(delta))) + 1;
				var weight = pt1*pt2*pt3;
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

	var normal_onesided_beta_imp = function(b0, b1, indiff, var_val, simulateResult, samples) {
		if (!samples) samples = 10000;
		var betas = [];
		var beta = 1; // precision/inverse-variance of the importance sampling distribution
		//var starttime = (new Date()).getTime();
		for (var i = 0;i < samples;i++) {
			var z = jStat.jStat.normal.sample(0,Math.sqrt(1/beta));
			z = Math.abs(z);

			var finished = false;
			var S_x = 0;
			var n = 0;
			var result = undefined;
			while (!finished) {
				n += 1;
				// pull xs from N(0,2*var_val)
				S_x += jStat.jStat.normal.sample(-z,Math.sqrt(2*var_val));

				var mle = S_x/n;
				// test on simplified boundaries
				if (mle < -indiff) {
					var L_na = 1;
				} else {
					var L_na = Math.exp( n*(S_x/n + indiff)*(S_x/n + indiff)/(4*var_val) );
				}
				if (L_na >= b0) {
					finished = true;
					result = "false"
				}
				if (mle > indiff) {
					var L_nb = 1;
				} else {
					var L_nb = Math.exp( n*(S_x/n - indiff)*(S_x/n - indiff)/(4*var_val) );
				}
				if (L_nb >= b1) {
					finished = true;
					result = "true";
				}
			}

			if (result == 'true') {
				var delta = (n + 2*var_val*beta)/(4*var_val);
				var pt1 = Math.sqrt(beta/(2*delta));
				var pt2 = Math.exp((-2*indiff*S_x + n*indiff*indiff)/(4*var_val) + (S_x*S_x)/(16*delta*var_val*var_val));
				var pt3 = 1 - jStat.jStat.erf(S_x/(4*var_val*Math.sqrt(delta)));
				var weight = pt1*pt2*pt3;
				betas.push(1/weight);
			} else {
				betas.push(0);
			}
		}
		//console.log("time:"+( (new Date()).getTime()-starttime ))
		//console.log("mean:"+mean(betas));
		//console.log("std_err:"+boot_std(betas,1000));
		return betas;
	}

	var normal_uv_onesided_LR_H0 = function(S_x, S_y, S_x2, S_y2, n, indiff) {
		if (n == 1) {
			// TODO : when n = 1, mle_llik is always -Infinity, so we have to enforce n > 1
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
			// TODO : when n = 1, mle_llik is always -Infinity, so we have to enforce n > 1
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

	var normal_kv_twosided_maxSamplesize = function(indiff, var_value, b0, b1) {
		var part = Math.sqrt(Math.log(b0)) + Math.sqrt(Math.log(b1));
		var max = part*part*(4*var_value)/(indiff*indiff);
		return max;
	}

	functions['normal_uv'] = {
		'two-sided' : {
			'l_an' : normal_uv_twosided_LR_H0,
			'l_bn' : normal_uv_twosided_LR_HA,
			//'alpha' : normal_twosided_alpha,
			'alpha' : normal_uv_twosided_alpha_imp3,
			//'beta' : normal_twosided_beta,
			'beta' : normal_uv_twosided_beta_imp2,
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
			'alpha' : normal_onesided_alpha_imp,
			'beta' : normal_onesided_beta_imp,
			'simulateH0' : normal_onesided_simulateH0,	
		}
	}
		  
	// precalculated thresholds
	thresholds['normal_uv'] = {
		'two-sided' : {
			0.05 : { // alpha
				0.10 : { // beta
					0.1 : { // indifference
						1 : [433, 8.85], // variance bound
					},
					0.05 : { // indifference
						1 : [482.5, 9.0], // variance bound, approximate
					},
					0.025 : { // indifference
						1 : [532, 9.1], // variance bound, needs to be tested more properly
					},
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
						1 : [137.4, 8.3] // variance
					},
					0.05 : { // indifference
						1 : [170.5, 8.55] // variance
					},
					0.025 : { // indifference
						1 : [205, 8.65] // variance
					},
					0.01 : {
						1 : [252.5, 8.6] // variance
					}
				}
			}
		},
		'one-sided' : {
			0.05 : { // alpha
				0.05 : { // beta
					0.1 : { // indifference
						1 : [51.3, 51.3], // variance
					},
					0.05 : { // indifference
						1 : [66, 66], // variance
					},
					0.025 : { // indifference
						1 : [81.3, 81.3], // variance
					},
					0.01 : { // indifference
						1 : [102.5, 102.5], // variance
					}
				}
			}
		}
	}

	tests['normal'] = normal_test;

	
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
			var lower_est_bc = optimize2d(lower_est, biasFun(), lower_est, 0.005, 16400, 590000, 0.02, undefined, undefined, false);
			var upper_est_bc = optimize2d(upper_est, biasFun(), upper_est, 0.005, 16400, 590000, 0.02, undefined, undefined, false);

			return [(lower_est_bc[0]-lower_est_bc[1]),(upper_est_bc[0]-upper_est_bc[1])];
		}

		// get estimate (only when test is done)
		  // use bias-reduction
		this.estimate = function() {
			if (!finished) {
				return undefined;
			}
			var ests = optimize2d([S_x/n_x, S_y/n_y], biasFun(), [S_x/n_x, S_y/n_y], 0.005, 16400, 590000, 0.02, undefined, undefined, false);
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
					var res = simulateResult([parameters[0], var_value], [parameters[1], var_value]);
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
	// debugging functions:
	  // test coverage of confidence intervals
}

if (typeof exports === 'object') {
	module.exports = new glr();
}
