options(dplyr.summarise.inform = FALSE)

`+` <- function(x, y) {
	if (typeof(x) == "character" |
			typeof(y) == "character")
		paste(as.character(x), as.character(y), sep = "")
	else
		.Primitive("+")(x, y)
}

ifelse = function(test, yes, no){
	if(test) return(yes)
	else return(no)
}

remove.na = function(df, col.names = colnames(df)){
	if(is.null(dim(df))) return(df[!is.na(df)])
	condition = T
	for(col.name in col.names){
		condition = condition & !is.na(df[[col.name]])
	}
	df[condition]
}

directory = function(..., sep = '/') {
	dirs = unlist(list(...))
	if(typeof(dirs) != 'character') stop('Directory names must be characters.')
	filepath = ''
	dirs = unlist(str_split(dirs, sep))
	for(dir in dirs){
		filepath = filepath + dir + .Platform$file.sep
		dir.create(filepath, showWarnings = F)
	}
	filepath
}

se <- function(v) {
	n.na = length(v[is.na(v)])
	n = length(v) - n.na
	if (n.na > 0) cat('sem() removed ' + n.na + ' NA elements (n = ' + n + ').\n')
	sd(v, na.rm = T) / sqrt(n)
}

ci = function(v, alpha = .95, null.as.0 = T) {
	m = mean(v, na.rm = T)
	if(length(v[!is.na(v)]) == 1){
		hw = 0
	}
	else{
		t = qt(1 - (1 - alpha) / 2, df = length(v[!is.na(v)]) - 1)
		hw = t * se(v)
		if(null.as.0 & is.null(hw)) hw = 0
	}
	list(
		'low' = m - hw,
		'high' = m + hw,
		'half_width' = hw,
		'hw' = hw
	)
}

signif.as.character = function(x, digits = 3){
	fun = function(x, digits){
		if (class(x) != 'numeric') return(x)
		if (x == 0) return('0.' + strrep('0', digits - 1))
		x.chr = as.character(signif(x, digits))
		significand = ''
		significant.digits = 0
		decimal = F
		exponential = F
		exponent = ''
		for(i in 1:nchar(x.chr)){
			letter = substring(x.chr, i, i)
			if(str_detect(letter, 'e')){
				exponential = T
				exponent = substring(x.chr, i, nchar(x.chr))
				break
			}
			if(str_detect(letter, '\\.')){
				decimal = T
			}
			significand = significand + letter
			if(str_detect(letter, ifelse(significant.digits == 0, '[1-9]', '[0-9]'))){
				significant.digits = significant.digits + 1
			}
		}
		if(exponential){
			if (significant.digits < digits){
				if(!decimal) significand = significand + '.'
				significand = significand + strrep('0', digits - significant.digits)
				x.chr = significand + exponent
			}
		}
		else if(significant.digits < digits){
			if(!decimal) x.chr = x.chr + '.'
			x.chr = x.chr + strrep('0', digits - significant.digits)
		}
		else if(significant.digits > digits){
			x.chr = substring(x.chr, 1, 1) + '.' + substring(x.chr, 2, digits) + 'e+' +
				ifelse(significant.digits - 1 < 10, '0', '') + (significant.digits - 1)
		}
		x.chr
	}
	sapply(x, fun, digits)
}

last.plot = function(gg = T){
	ifelse(gg, ggplot2::last_plot(), recordPlot())
}

save.plot = function(filepath = default.plot.path(), ..., plt = last.plot(), show = F){
	if(show) print(plt)
	ggsave(filepath, plt, ...)
	if(file.exists(filepath))	cat(filepath + ' saved.\n')
}

remove.outliters = function(dt, group.by, response, iterate = T){
	require(data.table)
	require(dplyr)
	summ = dt[, .(q1 = quantile(eval(parse(text = response)), .25),
								q3 = quantile(eval(parse(text = response)), .75),
								iqr = IQR(eval(parse(text = response)))),
						keyby = group.by]
	dt.rm = dt[summ, on = group.by][rfu > q1 - 1.5 * iqr & rfu < q3 + 1.5 * iqr] %>%
		mutate(q1 = NULL, q3 = NULL, iqr = NULL)
	if(iterate & nrow(dt) != nrow(dt.rm)) dt.rm = remove.outliters(dt.rm, group.by, response)
	dt.rm
}

lm.eq = function(model){
	a = unname(coef(model)[1])
	b = unname(coef(model)[2])
	r2 = summary(model)$adj.r.squared
	
	eq <- substitute(italic(y) == a~s~b~~italic(x)*","~~italic(R)^2~"="~r2, 
									 list(a = signif.as.character(a, digits = 4),
									 		 s = ifelse(b < 0, '-', '+'),
									 		 b = signif.as.character(b, digits = 4),
									 		 r2 = signif.as.character(r2, digits = 4)))
	as.character(as.expression(eq));
}

error.propagation = function(expr_text, ..., type = 'raw', alpha = 0.05){
	p = qpcR::propagate(parse(text = expr_text), cbind(...), type, alpha = alpha, plot = F)
	list(mean = p$summary$Prop[1],
			 sd = p$summary$Prop[2],
			 median = p$summary$Prop[3],
			 MAD = p$summary$Prop[4],
			 low = p$summary$Prop[5],
			 high = p$summary$Prop[6],
			 hw = p$summary$Prop[6] - p$summary$Prop[1])
}

#Convert Spark reading (RFU by emission wavelength in nm) into R data frame----
xl.to.dt = function(xl, mode = 'fluorescence'){
	xl = setDT(xl)
	is.wl = grepl(r'((\d+)nm)', colnames(xl))
	wls = as.integer(str_match(colnames(xl)[is.wl], '(\\d+)nm')[, 2])
	dt = xl[, !is.wl, with = F]
	result = NULL
	for(wl in wls) {
		dt.n = copy(dt)[, wl.nm := wl]
		if(mode %in% c('F', 'fluorescence')) dt.n[, rfu := xl[[wl + 'nm']]]
		if(mode %in% c('A', 'absorbance')) dt.n[, o.d. := xl[[wl + 'nm']]]
		result = rbind(result, dt.n, fill = T)
	}
	remove.na(result)
}

#Labels for salt solutes
solute.label = function(name){
	TeX(str_replace_all(str_replace_all(name, '(.+?)(\\d)', r'($\1_{\2}$)'), ' ', '~'))
}

#Wrap a list of plots into A4----
save.a4.wrap.plot = function(plots, filepath = waiver(), guides = 'collect', 
														 ncol = 3, nrow = NULL,
														 title = '', subtitle = '', ...){
	plt = wrap_plots(plots, guides = guides, ncol = ncol, nrow = nrow, ...) +
		plot_annotation(title = title,
										subtitle = subtitle)
	save.plot(filepath,	width = 210, height = 297, units = 'mm')
	plt
}

#Emission curve plots----
emission.plot = function(dt, group.by = c('solute', 'conc.M', 'wl.nm'),
												 mapping = aes(wl.nm, color = conc.M, fill = conc.M),
												 fill.label = 'Concentration [M]',
												 color.label = 'Concentration [M]', ...){
	if('conc.M' %in% names(dt)) dt = dt[order(-conc.M)]
	summ = dt[, .(m = mean(rfu), low = ci(rfu)$low, high = ci(rfu)$high), keyby = group.by]
	if(fill.label == 'Concentration [M]'){
		leg = formatC(signif(unique(dt$conc.M)), format = 'e', digits = 1)
		summ[, `:=`(conc.M = formatC(signif(conc.M), format = 'e', digits = 1))]
		if(length(unique(dt$conc.M)) == 1 & dt$conc.M[1] == 0) colors = 'gray'
		else colors = rainbow(length(unique(dt$conc.M)))
	}
	plt = summ %>% 
		ggplot(mapping) + 
		geom_vline(xintercept = 514, color = 'gray') +
		geom_line(aes(y = m)) +
		geom_ribbon(aes(ymin = low, ymax = high), alpha = 0.1, color = NA) +
		geom_text(x = 514.2, y = -Inf, hjust = 0, vjust = 0, label = '514', color = 'gray', size = 2) +
		labs(x = 'Emission wavelength [nm]', y = 'Fluorescence [RFU (95%CI)]', 
				 fill = fill.label, color = color.label, ...) 
	if(fill.label == 'Concentration [M]'){
		plt = plt + scale_color_manual(breaks = leg, values = colors) +
			scale_fill_manual(breaks = leg, values = colors)
	}
	plt
}

save.450.emission.plots = function(dt, category){
	lapply(unique(dt$solute), function(sol){
		lt = ifelse(sol == 'Water', 'longdash', 'solid')
		plt = emission.plot(dt[solute == sol], subtitle = solute.label(sol), linetype = lt)
		save.plot(directory('Figures/Ex 450nm emission spectrum/' + category) + sol + '.png')
		plt
	})
}


save.solution.ph.emission.plot = function(dt){
	emission.plot(dt, c('solution', 'ph', 'wl.nm'), 
								aes(wl.nm, color = solution, fill = solution, 
													 linetype = as.factor(ph)), 
								fill.label = '', color.label = '', linetype = 'pH',
								subtitle = 'Fluorescence by pH in AWE, PBS, Serum')
	save.plot(directory('Figures/PBS Serum AWE') + 'emission curve.png')
}

save.ph.emission.plot = function(dt){
	figs = list()
	n = length(unique(dt$ph))
	for(i in 1:n){
		p = unique(dt$ph)[i]
		figs[[i]] = emission.plot(dt[ph == p], c('ex.wl.nm', 'wl.nm'), 
															aes(wl.nm, color = as.factor(ex.wl.nm), fill = as.factor(ex.wl.nm)), 
															fill.label = TeX('$\\lambda_{ex}$'), 
															color.label = TeX('$\\lambda_{ex}$'),
															subtitle = 'ph = ' + p)
		save.plot(directory('Figures/By pH and excitation WL') + 'pH ' + p + '.png')
	}
	figs
}

#Absorption curve plots
absorption.plot = function(dt, control = NULL, group.by = c('solute', 'wl.nm', 'conc.M'), 
													 wl.range = c(370, 500),
													 mapping = aes(color = as.factor(conc.M), fill = as.factor(conc.M)),
													 fill.label = 'Concentration [M]',
													 color.label = 'Concentration [M]', ...){
	dt = copy(dt)[wl.nm >= wl.range[1] & wl.nm <= wl.range[2]]
	no.control = is.null(control)
	if(no.control){
		control = copy(dt)[, o.d. := 0]
		ax.y = 'Absorbance [OD (95%CI)]'
	}
	else{
		ax.y = TeX('Relative absorbance [OD - $OD_{0}$]') 
	}
	ctrl.summ = control[wl.nm >= wl.range[1] & wl.nm <= wl.range[2],
											.(ctrl = mean(o.d.), ctrl.ci = ci(o.d.)$hw), keyby = .(wl.nm)]
	dt = merge(dt, ctrl.summ, by = 'wl.nm')
	dt[, `:=`(o.d. = o.d. - ctrl)]
	if(fill.label == 'Concentration [M]'){
		leg = formatC(signif(unique(dt$conc.M)), format = 'e', digits = 1)
		dt = dt %>% mutate(conc.M = formatC(signif(conc.M), format = 'e', digits = 1))
		colors = rainbow(length(unique(dt$conc.M)))
	}
	if(!('solute' %in% names(dt))) dt[, solute := 'Water']
	summ = dt[, 
						.(m = mean(o.d.), 
							low = ifelse(is.null(control) | solute == 'Water', ci(o.d.)$low,
													 mean(o.d.) - hypot(ci(o.d.)$hw, sqrt(sum(ctrl.ci^2)))), 
							high = ifelse(is.null(control) | solute == 'Water', ci(o.d.)$high,
														mean(o.d.) + hypot(ci(o.d.)$hw, sqrt(sum(ctrl.ci^2))))), 
						by = c(group.by, 'solute')]
	plt = summ %>% 
		ggplot(mapping) + 
		geom_vline(xintercept = 405, color = 'red') +
		geom_vline(xintercept = 450, color = 'blue')
	if(!no.control) plt = plt + geom_hline(yintercept = 0)
	plt = plt +
		geom_line(aes(wl.nm, m)) +
		geom_ribbon(aes(x = wl.nm, ymin = low, ymax = high), alpha = 0.2, color = NA) +
		geom_text(x = 405.2, y = -Inf, hjust = 0, vjust = 0, label = '405', color = 'red', size = 2) +
		geom_text(x = 450.2, y = -Inf, hjust = 0, vjust = 0, label = '450', color = 'blue', size = 2) +
		
		labs(x = 'Wavelength [nm]', y = ax.y,
				 fill = fill.label, color = color.label, ...)
		if(fill.label == 'Concentration [M]'){
			plt = plt + scale_color_manual(breaks = leg, values = colors) +
				scale_fill_manual(breaks = leg, values = colors)
		}
		plt
}

save.absorption.plots = function(dt, control, category, wl.range = c(320, 1000)){
	for(sol in unique(dt$solute)){
		plt = absorption.plot(dt[solute == sol], control, wl.range = wl.range, 
													subtitle = solute.label(sol))
		save.plot(directory('Figures/Absorption/' + category) + sol + '.png')
	}
}

ph.absorption.plot = function(dt, control = NULL, wl.range = c(320, 550)){
	absorption.plot(dt, control, c('ph', 'wl.nm'), wl.range, 
									aes(color = as.factor(ph), fill = as.factor(ph)), 
									fill.label = 'pH', color.label = 'pH')
}

dilution.plot = function(dt, f0, colors, mode){
	is.individual = length(unique(dt$solute)) == 1
	y.label = switch(mode, 
									 'raw' = 'Fluorescence [RFU (95%CI)]',
									 'ratio' = TeX('Fluorescence [$\\frac{F-F_0}{F_0}$ (95%CI)]'),
									 'log' = TeX('Fluorescence [$\\log_2\\frac{F}{F_0}$ (95%CI)]'))
	conc = sort(unique(dt$conc.M), decreasing = F)
	summ = dt[, 
						.(m = mean(rfu), 
							low = switch(mode, 
													 'raw' = ci(rfu)$low,
													 'ratio' = error.propagation('(rfu - f0) / f0', rfu, f0)$loW,
													 'log' = mean(rfu) - sqrt(ci(rfu)$hw^2 + (ci(f0)$hw/mean(f0)/log(2))^2)),
							high = switch(mode, 
														'raw' = ci(rfu)$high,
														'ratio' = error.propagation('(rfu - f0) / f0', rfu, f0)$high,
														'log' = mean(rfu) + sqrt(ci(rfu)$hw^2 + (ci(f0)$hw/mean(f0)/log(2))^2))), 
						keyby = .(solute, conc.M)]
	plt = summ %>% ggplot(aes(x = conc.M, color = solute, fill = solute)) +
		geom_hline(aes(yintercept = 0), alpha = .5) +
		geom_point(data = dt, aes(y = rfu), alpha = .2, size = 1) +
		geom_line(aes(y = m), alpha = .5) +
		geom_errorbar(aes(ymin = low, ymax = high), width = 0, alpha = .5) +
		geom_ribbon(aes(ymin = low, ymax = high), alpha = 0.2, color = NA) +
		geom_vline(aes(xintercept = conc.awe[solute], color = solute), linetype = 'dashed', alpha = .5) +
		theme(axis.text.x = element_text(angle = 90)) +
		scale_x_continuous(breaks = conc, labels = formatC(signif(conc), format = 'e', digits = 1),
											 trans = 'log2') +
		scale_color_manual(breaks = unique(dt$solute), values = colors) +
		scale_fill_manual(breaks = unique(dt$solute), values = colors) +
		labs(x = 'Concentration [M]', y = y.label, color = '', fill = '')
	if(is.individual){
		plt = plt + labs(subtitle = solute.label(summ$solute[1])) + 
			guides(color = 'none', fill = 'none')
	}
	if(mode == 'raw'){
		plt = plt +	scale_y_continuous(trans = 'pseudo_log', breaks = log_breaks()(dt$rfu))
	}
	plt
}

save.dilution.plots = function(dt, category, ex.wl, control = NULL, mode = 'raw'){
	dir = directory('Figures/Concentration gradient @514nm/' +
										as.character(ex.wl) + 'nm/' + category + '/' + mode)
	dt = dt[order(solute)]
	solutes = unique(dt$solute)
	n = length(solutes)
	colors = setNames(rainbow(n), solutes)
	if(!is.null(control)){
		f0 = mean(control$rfu)
		if(mode == 'ratio'){
			dt$rfu = (dt$rfu - f0) / f0
		}
		else if(mode == 'log'){
			dt$rfu = log2(dt$rfu / f0)
		}
	}
	title = TeX(r'(Fluorescence at $\lambda_{em}$ = 514 nm ($\lambda_{ex}$ = )' + ex.wl + ' nm)')
	if(n > 1){
		dilution.plot(dt, control$rfu, colors, mode) + theme(legend.text = element_text(size = 6)) +
			labs(title = title, subtitle = category)
		save.plot(dir + 'Overview.png', width = 6.08, height = 4.24, units = 'in')
	}
	plots = lapply(unique(dt$solute), function(sol){
		fig = dilution.plot(dt[solute == sol], control$rfu, colors[sol], mode)
		save.plot(dir + directory('Individual') + sol + '.png')
		fig
	})
	save.a4.wrap.plot(plots, dir + 'Collective.png',
										title = title, subtitle = category, ncol = 3, nrow = 5)
}


find.peak.wl = function(dt, control){
	summ = dt[, .(m = mean(o.d.)), keyby = .(wl.nm)]
	c.summ = control[, .(m = mean(o.d.)), keyby = .(wl.nm)]
	summ$m = summ$m - c.summ$m
	cat('max: ' + summ$wl.nm[which.max(summ$m)] + '\n')
	cat('min: ' + summ$wl.nm[which.min(summ$m)] + '\n')
}

save.bar.plot = function(dt, conc, control, category, ex.wl){
	f0 = mean(control$rfu)
	dt$rfu = log2(dt$rfu / f0)
	solutes = unique(dt$solute)
	if(length(conc) == 1){
		dt = dt[conc.M == conc]
		filename = as.character(conc)
	}
	else if(length(conc) == length(solutes)){
		new.dt = NULL
		for(i in 1:length(solutes)){
			new.dt = rbind(new.dt,
										 dt[solute == solutes[i] & conc.M == conc[solutes[i]]])
		}
		dt = new.dt
		filename = as.character(min(conc)) + '—' + as.character(max(conc))
	}
	else{
		cat('Concentration length must be 1 or the same as the number of solutes.\n')
	}
	dt$solute = with(dt, reorder(solute, rfu, mean, decreasing = T))

	summ = dt[, .(low = ci(rfu)$low, high = ci(rfu)$high, rfu = mean(rfu)),
						keyby = .(solute, conc.M)][order(-rfu)] %>%
		.[, color := sapply(rfu, function(x) ifelse(x < 0, 'red', 'blue'))]
	just = sapply(summ$rfu, function(x){ifelse(x < 0, 0, 1)})
	plt = dt %>%
		ggbarplot('solute', 'rfu', fill = 'solute', palette = summ$color, orientation = 'horiz',
							add = c('mean_ci')) +
		geom_beeswarm(data = dt, mapping = aes(x = solute, y = rfu), size = 1, cex = 2, alpha = .33) +
		geom_text(data = summ,
							mapping = aes(y = 0, label = ' ' +
															formatC(signif(conc.M), format = 'e', digits = 1) + ' M '),
							hjust = just, vjust = 1, size = 2) +
		geom_hline(aes(yintercept = 0))	+
		labs(x = NULL, y = TeX(r'(Relative fluorescence $\log_2\frac{F}{F_0}$ (95%CI) )'),
				 title = 'Fluorescence (relative to water)',
				 subtitle = TeX(sprintf(r'($\lambda_{ex}=%d nm$, $\lambda_{em}=514 nm$)', ex.wl))) +
		scale_x_discrete(labels = parse(text = solute.label(summ$solute))) +
		theme(legend.position = "none")
	dir = directory('Figures/Concentration gradient @514nm (Normalized, log2)/Summary/' +
										as.character(ex.wl) + 'nm/' + category)
	save.plot(dir + filename +' M.png')
}

save.solution.ph.heatmap = function(dt){
	dt[wl.nm == 514, .(m = mean(rfu)), keyby = .(solution, ph)] %>%
		ggplot(aes(x = solution, y = as.factor(ph), fill = m)) + 
		geom_tile() +
		geom_text(aes(label = round(m)), color = 'white') + 
		labs(x = NULL, y = 'pH', fill = 'RFU', subtitle = 'Fluorescence [RFU] @514 nm')
	save.plot(directory('Figures/PBS Serum AWE') + 'heatmap.png')
}


closest.conc = function(v.conc, v.dilutions){
	v = sapply(v.conc, function(x){
		v.dilutions[which.min(abs(v.dilutions - x))]
	})
	names(v) = names(v.conc)
	v
}
