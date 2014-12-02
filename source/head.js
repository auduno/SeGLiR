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
	var optimize2d = function(vec, fun, init_points, epsilon, nsamples, max_samples, gradient_evaluation_width, lower_limit, upper_limit, verbose, debug) {
		if (typeof(verbose) === 'undefined') {
			verbose = true;
		}
		if (typeof(debug) === 'undefined') {
			debug = false;
		}

		// clone variables
		vec = [vec[0],vec[1]];
		init_points = [init_points[0],init_points[1]];

		var gradient, point, next_point;
		var samples = 0;
		var diff = Infinity;
		var est_point = init_points;
		if (verbose) console.log("Initial estimate : "+est_point);
		var iteration = 0;
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
				if (debug) console.log("alpha : "+curpoints[0]+" +- "+(4*std(curval[0])/Math.sqrt(curval[0].length)));
				if (debug) console.log("beta : "+curpoints[1]+" +- "+(4*std(curval[1])/Math.sqrt(curval[1].length)));

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

				if (debug) console.log("checked current estimate, samples:"+samples)
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
				
				if (debug) console.log("length samples : "+l_samples[0].length);

				samples += num_samples;
				if (debug) console.log(samples);

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

				if (debug) console.log("b1p1 gradient : "+b1p1_gradient_mean+" +- "+(4*Math.sqrt(b1p1_var)) );
				if (debug) console.log("b1p2 gradient : "+b1p2_gradient_mean+" +- "+(4*Math.sqrt(b1p2_var)) );
				if (debug) console.log("b2p1 gradient : "+b2p1_gradient_mean+" +- "+(4*Math.sqrt(b2p1_var)) );
				if (debug) console.log("b2p2 gradient : "+b2p2_gradient_mean+" +- "+(4*Math.sqrt(b2p2_var)) );
				
				enoughSamples = true;
				if (sign(b1p1_gradient_mean+4*Math.sqrt(b1p1_var)) != sign(b1p1_gradient_mean-4*Math.sqrt(b1p1_var))) enoughSamples = false;
				if (sign(b2p2_gradient_mean+4*Math.sqrt(b2p2_var)) != sign(b2p2_gradient_mean-4*Math.sqrt(b2p2_var))) enoughSamples = false;
				//if (sign(b1p2_gradient_mean+4*Math.sqrt(b1p2_var)) != sign(b1p2_gradient_mean-4*Math.sqrt(b1p2_var))) enoughSamples = false;
				//if (sign(b2p1_gradient_mean+4*Math.sqrt(b2p1_var)) != sign(b2p1_gradient_mean-4*Math.sqrt(b2p1_var))) enoughSamples = false;
				
				if (debug) console.log("b1p1*b2p2:"+(b1p1_gradient_mean*b2p2_gradient_mean));
				if (debug) console.log("b1p2*b2p1:"+(b1p2_gradient_mean*b2p1_gradient_mean));
				if (debug) console.log("b2p2*(v0-c0):"+( b2p2_gradient_mean*(vec[0]-curpoints[0]) ));
				if (debug) console.log("b1p2*(v1-c1):"+( b1p2_gradient_mean*(vec[1]-curpoints[1]) ));

				i += 1;
				if (debug) console.log("getting gradients, samples:"+samples)
			}

			// extrapolate where point lies with simple linear function
			var mult = 1/(b1p1_gradient_mean*b2p2_gradient_mean - b1p2_gradient_mean*b2p1_gradient_mean);
			next_point = [mult,mult];
			next_point[0] *= ( b2p2_gradient_mean*(vec[0]-curpoints[0]) - b1p2_gradient_mean*(vec[1]-curpoints[1]) );
			next_point[1] *= ( b1p1_gradient_mean*(vec[1]-curpoints[1]) - b2p1_gradient_mean*(vec[0]-curpoints[0]) );

			// calculate difference between new point and estimated point
			//diff = Math.sqrt(next_point[0]*next_point[0] + next_point[1]*next_point[1]);

			//next_point = est_point + 0.5*(next_point-est_point);
			if (debug) console.log(next_point);
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

			iteration += 1;
			if (verbose) console.log("Iteration "+iteration+" estimate : "+est_point);
		}

		if (verbose && samples > max_samples) {
			console.log("Stopped estimation due to sample limit reached. Estimate did not converge.")
		}

		if (verbose) console.log("Final estimate : "+est_point);
		
		// calculate estimate of final value
		var curval = fun(est_point, nsamples);
		curpoints = [mean(curval[0]), mean(curval[1])];
		if (debug) console.log("alpha : "+curpoints[0]+" +- "+(4*std(curval[0])/Math.sqrt(nsamples)));
		if (debug) console.log("beta : "+curpoints[1]+" +- "+(4*std(curval[1])/Math.sqrt(nsamples)));

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
