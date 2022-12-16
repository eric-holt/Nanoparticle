packages = c('data.table', 'ggbeeswarm', 'ggh4x', "ggpubr", 'ggpmisc', 'ggrepel', 'latex2exp', 'NCmisc', "patchwork", 'pracma', 'propagate', 'qpcR', 'rlang', 'rstatix', 'scales', "tidyverse", 'readxl')
for (package in packages) {
	if (! package %in% installed.packages())
		install.packages(package, dependencies = T)
	require(package, character.only = T)
}

source('Functions.R')

theme_replace(plot.title = element_text(hjust = .5, size = 10),
							plot.subtitle = element_text(hjust = .5, size = 8),
							plot.caption = element_text(hjust = 1),
							plot.tag = element_text(face = 'italic'),
							legend.title = element_text(size = 7),
							legend.text = element_text(size = 7), 
							axis.title = element_text(size = 8),
							axis.text = element_text(size = 7))

#Read from Excel----
hpts.ex.em.xl = read_excel('Data/R/HPTS & CHITOSAN+HPTS Ex Em.xlsx', sheet = 'HPTS DYE')
chit.ex.em.xl = read_excel('Data/R/HPTS & CHITOSAN+HPTS Ex Em.xlsx', sheet = 'CHITOSAN + HPTS')
aa.o.450.xl = read_excel("Data/R/AA Organic.xlsx", sheet = "Fluorescence Ex 450nm")
salt.450.xl = read_excel("Data/R/Salts.xlsx", sheet = "Fluorescence Ex 450nm")
aa.o.405.xl = read_excel("Data/R/AA Organic.xlsx", sheet = "Fluorescence Ex 405nm")
salt.405.xl = read_excel("Data/R/Salts.xlsx", sheet = "Fluorescence Ex 405nm")
aa.o.abs.xl = read_excel("Data/R/AA Organic.xlsx", sheet = "Absorption")
salt.abs.xl = read_excel("Data/R/Salts.xlsx", sheet = "Absorption")
ph.emission.xl = read_excel("Data/R/pH PBS.xlsx", sheet = 'Emission')
ph.absorpnion.xl = read_excel('Data/R/pH PBS.xlsx', sheet = 'Absorption')
ph.emission.450.xl = read_excel("Data/R/pH PBS.xlsx", sheet = 'HPTS CHIT dilution 25ul + 50ul ')
ph.pbs.serum.awe.xl = read_excel('Data/R/pH PBS Serum AWE.xlsx', sheet = 'Fluorescence PBS Serum AWE')
aa.conc.xl = read_excel('Data/R/AWE Concentrations.xlsx', sheet = 'Amino acids')
salt.conc.xl = read_excel('Data/R/AWE Concentrations.xlsx', sheet = 'Salts')
o.conc.xl = read_excel('Data/R/AWE Concentrations.xlsx', sheet = 'Other')
hpts.chit.conc.xl = read_excel('Data/R/HPTS conc in Chitosan ADAM.xlsx', 
															 sheet = 'HPTS conc in Chitosan')
hpts.chit.ex.ph.xl = read_excel('Data/R/Excitation Emission HPTS-CHIT 10 nm steps.xlsx', 
																sheet = 'HPTS-CHIT')

#data tables
#pH----
ph.emission = xl.to.dt(ph.emission.xl)
ph.emission.450 = xl.to.dt(ph.emission.450.xl)
ph.absorpnion = xl.to.dt(ph.absorpnion.xl, 'A')
ph.pbs.serum.awe = xl.to.dt(ph.pbs.serum.awe.xl)

#Amino acids, salts, organic substances---
aa.o.450 = xl.to.dt(aa.o.450.xl)
salt.450 = xl.to.dt(salt.450.xl)
aa.o.405 = xl.to.dt(aa.o.405.xl)
salt.405 = xl.to.dt(salt.405.xl)
aa.o.abs = xl.to.dt(aa.o.abs.xl, 'A')
salt.abs = xl.to.dt(salt.abs.xl, 'A')

#HPTS, chitosan----
chit.ex.em = xl.to.dt(chit.ex.em.xl)
hpts.ex.em = xl.to.dt(hpts.ex.em.xl)
hpts.chit.ex.em = remove.na(rbind(chit.ex.em %>% mutate(solute = 'HPTS-Chitosan'),
																	hpts.ex.em %>% mutate(solute = 'HPTS')))
hpts.chit.conc = xl.to.dt(hpts.chit.conc.xl)
hpts.chit.ex.ph = xl.to.dt(hpts.chit.ex.ph.xl) #<from a different data set

#ex. 450nm emission data----
aa.450 = aa.o.450[nchar(solute) == 3]
o.450 = aa.o.450[nchar(solute) > 3]
salt.450.control = salt.450[grepl('no_HPTS', solute) | solute == 'Water']
salt.450.data = salt.450[!grepl('no_HPTS', solute) & solute != 'Water']

#ex. 450nm concentration data at em. 514nm----
aa.450.514 = aa.450[wl.nm == 514]
o.450.514 = o.450[wl.nm == 514 & solute != 'Water' & solute != 'Albumin']
alb.450.514 = o.450[wl.nm == 514 & solute == 'Albumin']
salt.450.514 = salt.450.data[wl.nm == 514]

#ex. 405nm concentration data at em. 514nm----
aa.405.514 = aa.o.405[nchar(solute) == 3]
o.405.514 = aa.o.405[nchar(solute) > 3 & solute != 'Water' & solute != 'Albumin']
alb.405.514 = aa.o.405[solute == 'Albumin']
salt.405.514 = salt.405[!grepl('no_HPTS', solute) & solute != 'Water']

#ex. 450nm, em. 514nm, Water----
aa.o.water.450.514 = aa.o.450[wl.nm == 514 & solute == 'Water']
salt.water.450.514 = aa.o.405[wl.nm == 514 & solute == 'Water']

#ex. 405nm, em. 514nm, Water----
aa.o.water.405.514 = aa.o.405[solute == 'Water']
salt.water.405.514 = salt.405[solute == 'Water']

#absorption----
aa.o.abs.water = aa.o.abs[solute == 'Water']
aa.abs = aa.o.abs[nchar(solute) == 3]
o.abs = aa.o.abs[nchar(solute) > 3 & solute != 'Water']
salt.abs.water = salt.abs[solute == 'Water']
salt.abs.control = salt.abs[grepl('no_HPTS', solute)][order(solute)]
salt.abs.data = salt.abs[!grepl('no_HPTS', solute) & solute != 'Water'][order(solute)]#HPTS vs HPTS-chitosan plot----
#Emission plots
p1 = hpts.chit.ex.em[ex.wl.nm <= 480 & wl.nm >= 460 & wl.nm <= 580] %>%
	emission.plot(c('ex.wl.nm', 'wl.nm', 'solute'), aes(wl.nm, color = as.factor(ex.wl.nm)),
								color.label = TeX('$\\lambda_{ex}$'), fill.label = TeX('$\\lambda_{ex}$')) +
	facet_nested(. ~ solute, scales = 'free', independent = 'y') +
	labs(y = 'Fluorescence [RFU]') +
	theme(legend.text = element_text(size = 6))

#Peak intensity
p2 = hpts.chit.ex.em[wl.nm == 514, 
										 .(m = mean(rfu), low = ci(rfu)$low, high = ci(rfu)$high)
										 , keyby = .(ex.wl.nm, solute)] %>% 
	ggplot(aes(ex.wl.nm, m, color = solute, fill = solute)) + 
	geom_vline(xintercept = 405, color = 'red') +
	geom_vline(xintercept = 450, color = 'blue') +
	geom_text(x = 405.2, y = -Inf, hjust = 0, vjust = 0, label = '405', color = 'red', size = 2) +
	geom_text(x = 450.2, y = -Inf, hjust = 0, vjust = 0, label = '450', color = 'blue', size = 2) +
	geom_line() + geom_point() +
	geom_ribbon(aes(ymin = low, ymax =high), alpha = 0.1, color = NA) +
	labs(x = 'Excitation wavelength [nm]', y = 'Fluorescence [RFU]', 
			 color = '', fill = '',
			 subtitle = TeX(r'(Fluorescence at $\lambda_{em}$ = 514 nm)')) +
	scale_x_continuous(breaks = unique(hpts.chit.ex.em$ex.wl.nm))

p1 / p2 + plot_annotation(tag_levels = 'A')
save.plot(directory('Figures/Report') + 'HPTS vs HPTS-chitosan.png', 
					width = 6.31, height = 5, units = 'in')

#pH PBS----
#Emission spectra
# dt = rbind(ph.emission, hpts.chit.ex.ph)
dt = ph.emission
# dt$ph = factor(dt$ph, levels = c(5, 9, 5.5, 9.5))
p1 = emission.plot(dt[wl.nm >= 450 & wl.nm <= 580], c('ph', 'ex.wl.nm', 'wl.nm'), 
									 aes(wl.nm, color = as.factor(ex.wl.nm), fill = as.factor(ex.wl.nm)),
									 fill.label = TeX('$\\lambda_{ex}$'), color.label = TeX('$\\lambda_{ex}$')) +
	facet_wrap(~ paste0('pH ', ph), nrow = 1) +
	guides(color = guide_legend(ncol = 2), fill = guide_legend(ncol = 2))

#Peak line graph with the pH 5 + 9 data
p2 = ph.emission[wl.nm == 514 | wl.nm == 512, 
								 .(m = mean(rfu), low = ci(rfu)$low, high = ci(rfu)$high), 
								 by = .(ex.wl.nm, ph, wl.nm)] %>%
	ggplot(aes(ex.wl.nm, m, color = as.factor(ph), fill = as.factor(ph))) + 
	geom_vline(xintercept = 405, color = 'red') +
	geom_vline(xintercept = 450, color = 'blue') +
	geom_text(x = 405.2, y = -Inf, hjust = 0, vjust = 0, label = '405', color = 'red', size = 2) +
	geom_text(x = 450.2, y = -Inf, hjust = 0, vjust = 0, label = '450', color = 'blue', size = 2) +
	geom_line() +  geom_point(size = .5) +
	geom_ribbon(aes(ymin = low, ymax =high), alpha = 0.2, color = NA) +
	labs(x = 'Excitation wavelength [nm]', y = 'Fluorescence [RFU (95%CI)]', 
			 color = 'pH', fill = 'pH',
			 subtitle = TeX(r'($\lambda_{em}$ = 514 nm)')) +
	scale_x_continuous(breaks = unique(ph.emission$ex.wl.nm)) +
	theme(axis.text.x = element_text(angle = 90))

#Peak line graph with the pH 5.5 + 9.5 data
p3 = hpts.chit.ex.ph[wl.nm == 514 | wl.nm == 512, 
										 .(m = mean(rfu), low = ci(rfu)$low, high = ci(rfu)$high), 
										 by = .(ex.wl.nm, ph, wl.nm)] %>%
	ggplot(aes(ex.wl.nm, m, color = as.factor(ph), fill = as.factor(ph))) + 
	geom_vline(xintercept = 405, color = 'red') +
	geom_vline(xintercept = 450, color = 'blue') +
	geom_text(x = 405.2, y = -Inf, hjust = 0, vjust = 0, label = '405', color = 'red', size = 2) +
	geom_text(x = 450.2, y = -Inf, hjust = 0, vjust = 0, label = '450', color = 'blue', size = 2) +
	geom_line() + geom_point(size = .5) +
	geom_ribbon(aes(ymin = low, ymax =high), alpha = 0.2, color = NA) +
	labs(x = 'Excitation wavelength [nm]', y = 'Fluorescence [RFU]', 
			 color = 'pH', fill = 'pH',
			 subtitle = TeX(r'($\lambda_{em}$ = 514 nm)')) +
	scale_x_continuous(breaks = unique(hpts.chit.ex.ph$ex.wl.nm)) +
	theme(axis.text.x = element_text(angle = 90))

#Plots with pH 5–9.5 at ex.450
#Emission
p4 = emission.plot(ph.emission.450[wl.nm <= 560], c('ph', 'wl.nm'), 
									 aes(wl.nm, color = as.factor(ph), fill = as.factor(ph)), 
									 fill.label = 'pH', color.label = 'pH') + 
	labs(subtitle = TeX(r'($\lambda_{ex} = 450$ nm)')) +
	guides(color = guide_legend(ncol = 2), fill = guide_legend(ncol = 2))

#Peak transition over pH
p5 = ph.emission.450[wl.nm == 514, .(rfu = mean(rfu), low = ci(rfu)$low, high = ci(rfu)$high), 
										 keyby = c('ph', 'wl.nm')] %>%
	ggplot(aes(ph, rfu)) +
	geom_point(size = .5) + 
	geom_errorbar(aes(ymin = low, ymax = high), alpha = .5, width = .1) + 
	labs(x = 'pH', y = 'Fluorescence [RFU (95%CI)]', 
			 subtitle = TeX(r'($\lambda_{ex} = 450$ nm, $\lambda_{em} = 514$ nm)'))

#Absorption spectra
p6 = 	absorption.plot(ph.absorpnion, group.by = c('ph', 'wl.nm'), 
											mapping = aes(color = as.factor(ph), fill = as.factor(ph)), 
											fill.label = 'pH', color.label = 'pH') 
	

#Peak transition over pH
p7 = ph.absorpnion[wl.nm == 450, 
									 .(o.d. = mean(o.d.), low = ci(o.d.)$low, high = ci(o.d.)$high),
									 keyby = .(ph, wl.nm)] %>%
	ggplot(aes(ph, o.d.)) +
	geom_point(size = .5) + 
	geom_errorbar(aes(ymin = low, ymax = high), alpha = .5, width = .1) + 
	labs(x = 'pH', y = 'Absorbance [OD (95%CI)]', 
			 subtitle = TeX(r'($\lambda_{abs} = 450$ nm)'))

(p5 + p4 + theme(axis.title.y = element_blank())) / 
	(p7 + p6 + theme(axis.title.y = element_blank())) / 
	(p2 + p1 + theme(axis.title.y = element_blank())) + plot_annotation(tag_levels = 'A')
save.plot(directory('Figures/Report') + 'pH.png', 
					width = 6.31, height = 5.5, units = 'in')


#Concentration comparison plot----
pbs.rfu = hpts.chit.conc[wl.nm == 514 & solute == 'PBS', .(mean = mean(rfu), ci = ci(rfu)$hw)]
hpts.chit.conc.summ = 
	hpts.chit.conc[wl.nm == 514 & solute != 'PBS'] %>%
	mutate(rfu = rfu - pbs.rfu$mean) %>%
	.[, .(m = mean(rfu), 
				low = mean(rfu) - hypot(ci(rfu)$hw, pbs.rfu$ci), 
				high = mean(rfu) + hypot(ci(rfu)$hw, pbs.rfu$ci),
				hw = hypot(ci(rfu)$hw, pbs.rfu$ci)), 
		keyby = .(solute, conc.ppm)]
hpts = hpts.chit.conc.summ[solute == 'HPTS']
chit = hpts.chit.conc.summ[solute == 'HPTS-Chitosan']
#Linear regression + fitting
# dt.reg: 12 HPTS data points used for regression
dt.reg = hpts.chit.conc[wl.nm == 514 & solute == 'HPTS' & conc.ppm > .4 & conc.ppm < 4]
# dt.chit: 12 HPTS-chitosan data points used for fitting
dt.chit = hpts.chit.conc[wl.nm == 514 & solute == 'HPTS-Chitosan']
lm.reg = lm(rfu ~ conc.ppm, dt.reg)
dt.fit.x = data.table(conc.ppm = dt.chit$conc.ppm, 
											fit = sapply(dt.chit$rfu, function(y){
												(y - unname(coef(lm.reg)[1])) / unname(coef(lm.reg)[2])
											}))
dt.fit.x[, hpts.percent := fit / conc.ppm * 100]
pred = predict(lm.reg, newdata = data.table(conc.ppm = dt.fit.x$fit), interval = 'prediction')
dt.fit.x[, ci := sapply(pred[,3], function(y){
	(y - unname(coef(lm.reg)[1])) / unname(coef(lm.reg)[2])
}) / dt.chit$conc.ppm * 100 - hpts.percent]
dt.fit = dt.fit.x[, .(x = mean(fit),
											hpts.percent = mean(hpts.percent),
											ci = sqrt(sum(ci^2))), keyby = .(conc.ppm)]

mean.percent = mean(dt.fit$hpts.percent)
propagated.ci.hw = sqrt(sum(dt.fit$ci^2))
#Plot
plt1 = hpts.chit.conc.summ %>%
	ggplot(aes(conc.ppm, m, color = solute, fill = solute)) +
	geom_line() +
	geom_point() +
	geom_ribbon(aes(ymin = low, ymax = high), alpha = .2, color = NA) +
	labs(x = TeX('Concentration \\[ppm\\] ($log_{2}$)'), 
			 y = TeX('Fluorescence \\[RFU 95%CI\\] ($log_{2}$)'),
			 color = '', fill = '',
			 subtitle = TeX(r'(HPTS content in HPTS-chitosan based on fluorescence ($\lambda_{ex}$ = 450 nm, \lambda_{em}$ = 514 nm))')) +
	scale_x_continuous(breaks = unique(hpts.chit.conc$conc.ppm), 
										 labels = signif.as.character(unique(hpts.chit.conc$conc.ppm)),
										 trans = 'log2') +
	scale_y_continuous(
		breaks = unique(chit$m),
		labels = signif.as.character(unique(chit$m), 2),
		trans = 'log2')
for(i in 1:4){
	plt1 = plt1 + 	
		geom_text(x = log2(dt.fit$x[i]), y = log2(0), hjust = 0, vjust = 0.5, 
							angle = 90, color = 'purple', size = 2.5,
							label = signif.as.character(dt.fit$x[i]) + '\n(' + 
								format(round(dt.fit$hpts.percent[i], 2), nsmall = 2) + '±' + 
								format(round(dt.fit$ci[i], 2), nsmall = 2) + '%)')
}
plt1 = plt1 +
	geom_segment(aes(x = dt.fit$x[1], y = chit$m[1], xend = chit$conc.ppm[1], yend = chit$m[1]),
							 linetype = 'dotted') +
	geom_segment(aes(x = dt.fit$x[2], y = chit$m[2], xend = chit$conc.ppm[2], yend = chit$m[2]),
							 linetype = 'dotted') +
	geom_segment(aes(x = dt.fit$x[3], y = chit$m[3], xend = chit$conc.ppm[3], yend = chit$m[3]),
							 linetype = 'dotted') +
	geom_segment(aes(x = dt.fit$x[4], y = chit$m[4], xend = chit$conc.ppm[4], yend = chit$m[4]),
							 linetype = 'dotted') +
	geom_segment(aes(x = dt.fit$x[1], y = 0, xend = dt.fit$x[1], yend = chit$m[1]), 
							 linetype = 'dotted') +
	geom_segment(aes(x = dt.fit$x[2], y = 0, xend = dt.fit$x[2], yend = chit$m[2]), 
							 linetype = 'dotted') +
	geom_segment(aes(x = dt.fit$x[3], y = 0, xend = dt.fit$x[3], yend = chit$m[3]), 
							 linetype = 'dotted') +
	geom_segment(aes(x = dt.fit$x[4], y = 0, xend = dt.fit$x[4], yend = chit$m[4]), 
							 linetype = 'dotted')
plt1 = plt1 +
	geom_text(x = log2(min(chit$conc.ppm)), y = log2(min(hpts$m)), hjust = 0, vjust = 0, 
						color = 'purple', size = 2.5,
						label = 'Mean: ' + format(round(mean.percent, 2), nsmall = 2) + '±' + 
							format(round(propagated.ci.hw, 2), nsmall = 2) + '% HPTS') +
	geom_text(x = log2(0.3), y = log2(15000), hjust = 0, vjust = 0, 
						label = lm.eq(lm.reg), parse = T, 
						color = 'blue', alpha = .5, size = 2.5) +
	geom_smooth(aes(y = rfu, color = NULL, fill = NULL), data = dt.reg, method = 'lm', 
							color = 'blue', fill = 'blue', linewidth = .2)
save.plot(directory('Figures/Report') + 'HPTS in chitosan composite.png', 
					width = 6.31, height = 4, units = 'in')


#Compounds----
#450 Emission (AA organic)
dt = copy(aa.o.450)
dt$solute = factor(dt$solute,
									 levels = c(unique(aa.450$solute), unique(o.450$solute)))
emission.plot(dt[wl.nm <= 550]) + facet_wrap2(~ solute, scales = 'free', ncol = 4)
save.plot(directory('Figures/Report') + '450 Emission (AA + organic).png', 
					width = 6.31, height = 5, units = 'in')

#450 Emission (salts)
dt = copy(salt.450)
dt$solute = 
	factor(dt$solute,
				 levels = c(unique(salt.450.data$solute), 
				 					 unique(salt.450.control[solute != 'Water']$solute), 'Water'),
				 labels = c(solute.label(unique(salt.450.data$solute)),
				 					 solute.label(str_replace(unique(salt.450.control[solute != 'Water']$solute),
				 					 												 '_no_HPTS', '_ctrl')), 'Water'))
emission.plot(dt[wl.nm <= 550]) + facet_wrap(~ solute, scales = 'free', ncol = 4, 
																						 labeller = label_parsed)
save.plot(directory('Figures/Report') + '450 Emission (Salts).png', 
					width = 6.31, height = 7, units = 'in')

#Ex 450 Em 514 dilution
dilution.plot.report = function(dt, control){
	f0 = control$rfu
	dt$solute.label = factor(dt$solute, 
													 levels = unique(dt$solute), 
													 labels = solute.label(unique(dt$solute)))
	dt$rfu = log2(dt$rfu/mean(f0))
	summ = dt[, 
						.(rfu = mean(rfu), 
							low = mean(rfu) - sqrt(ci(rfu)$hw^2 + (ci(f0)$hw/mean(f0)/log(2))^2),
							high = mean(rfu) + sqrt(ci(rfu)$hw^2 + (ci(f0)$hw/mean(f0)/log(2))^2)), 
						keyby = .(solute, solute.label, conc.M)]
	summ.closest = summ[conc.M == conc[solute]]
	summ %>% ggplot(aes(conc.M, rfu)) +
		geom_hline(aes(yintercept = 0), alpha = .5) +
		geom_point(data = dt, alpha = .2, size = .5) +
		geom_line(alpha = .5) +
		geom_errorbar(aes(ymin = low, ymax = high), width = 0, alpha = .5) +
		geom_ribbon(aes(ymin = low, ymax = high), alpha = 0.2, color = NA) +
		geom_vline(aes(xintercept = conc.awe[solute]), linetype = 'dashed', color = 'red', alpha = .5) +
		geom_point(data = summ.closest, color = 'red', size = .5) +
		geom_errorbar(aes(ymin = low, ymax = high), data = summ.closest, width = 0, color = 'red') +
		facet_wrap(~ solute.label, ncol = 5, scales = 'free_y', labeller = label_parsed) +
		theme(axis.text.x = element_text(angle = 90)) +
		scale_x_continuous(breaks = unique(dt$conc.M)[c(T, F)], 
											 labels = formatC(signif(unique(dt$conc.M)[c(T, F)]), format = 'e', digits = 1),
											 trans = 'log2') +
		labs(x = 'Concentration [M]', y = TeX('Fluorescence [$\\log_2\\frac{F}{F_0}$ (95%CI)]'))
}
p1 = dilution.plot.report(aa.450.514, aa.o.water.450.514)
p2 = dilution.plot.report(o.450.514, aa.o.water.450.514)
p3 = dilution.plot.report(alb.450.514, aa.o.water.450.514)
p4 = dilution.plot.report(salt.450.514, salt.water.450.514)
p1$labels$x = p2$labels$x = p3$labels$x = '' 
p1$labels$y = p3$labels$y = p4$labels$y = ''
p1 / (p2 + p3 + plot_layout(widths = c(6, 1))) / p4 + plot_layout(heights = c(4, 1, 4))
save.plot(directory('Figures/Report') + '514 nm Dilution (450 nm).png', 
					width = 6.31, height = 7.5, units = 'in')

#Ex 405 Em 514 dilution
p1 = dilution.plot.report(aa.405.514, aa.o.water.405.514)
p2 = dilution.plot.report(o.405.514, aa.o.water.405.514)
p3 = dilution.plot.report(alb.405.514, aa.o.water.405.514)
p4 = dilution.plot.report(salt.405.514, salt.water.405.514)
p1$labels$x = p2$labels$x = p3$labels$x = '' 
p1$labels$y = p3$labels$y = p4$labels$y = ''
p1 / (p2 + p3 + plot_layout(widths = c(6, 1))) / p4 + plot_layout(heights = c(4, 1, 4))
save.plot(directory('Figures/Report') + '514 nm Dilution (405 nm).png', 
					width = 6.31, height = 7.5, units = 'in')

#AA organic absorption
dt = copy(aa.o.abs)
dt$solute = factor(dt$solute,
									 levels = c(unique(aa.450$solute), unique(o.450$solute)),
									 labels = solute.label(c(unique(aa.450$solute), unique(o.450$solute))))
absorption.plot(dt, aa.o.abs.water) +
	facet_wrap(~solute, scales = 'free_y', ncol = 3, labeller = label_parsed)
save.plot(directory('Figures/Report/') + 'Absorption (AA + organic).png',
					width = 6.31, height = 7, units = 'in')

#Salt absorption
dt = copy(salt.abs)
dt$solute = 
	factor(dt$solute,
				 levels = c(unique(salt.abs.data$solute), unique(salt.abs.control$solute), 'Water'),
				 labels = 
				 	solute.label(c(unique(salt.abs.data$solute), 
				 								 str_replace(unique(salt.abs.control$solute), 'no_HPTS', 'ctrl'), 'Water')))
absorption.plot(dt, salt.abs.water) +
	facet_wrap(~solute, scales = 'free_y', ncol = 4, labeller = label_parsed)
save.plot(directory('Figures/Report/') + 'Absorption (Salts).png',
					width = 6.31, height = 7, units = 'in')

#Bar plots
dt.bar.summary = function(dt, control, conc){
	dt = copy(dt)
	f0 = control$rfu
	dt$solute = factor(dt$solute, 
										 levels = unique(dt$solute))
	dt$rfu = log2(dt$rfu / mean(f0))
	solutes = unique(dt$solute)
	if(length(conc) == length(solutes)){
		new.dt = NULL
		for(i in 1:length(solutes)){
			new.dt = rbind(new.dt,
										 dt[solute == solutes[i] & conc.M == conc[solutes[i]]] %>% 
										 	 	.[, `:=`(ctrl = mean(f0), ctrl.ci = ci(f0)$hw)])
		}
		dt = new.dt
	}
	else{
		cat('Concentration length must be 1 or the same as the number of solutes.\n')
	}
	dt
}

bar.plot.report = function(dt){
	dt = copy(dt)
	dt$solute = with(dt, reorder(solute, rfu, mean, decreasing = T))
	summ = dt[,
						.(rfu = mean(rfu),
							low = mean(rfu) - sqrt(ci(rfu)$hw^2 + (mean(ctrl.ci)/mean(ctrl)/log(2))^2),
							high = mean(rfu) + sqrt(ci(rfu)$hw^2 + (mean(ctrl.ci)/mean(ctrl)/log(2))^2)),
						by = .(solute, conc.M)] %>%
		.[, color := sapply(rfu, function(x) ifelse(x < 0, '1', '2'))] %>%
		.[order(-rfu)]
	just = sapply(summ$rfu, function(x){ifelse(x < 0, 0, 1)})
	summ %>% ggplot(aes(solute, rfu)) +
		geom_hline(yintercept = 0, alpha = .5) +
		geom_col(aes(fill = color)) +
		geom_errorbar(aes(ymin = low, ymax = high), width = .2, alpha = .5) +
		geom_beeswarm(data = dt, size = 1, alpha = .5) +
		geom_text(mapping = aes(y = 0, angle = 90, label = ' ' +
															formatC(signif(conc.M), format = 'e', digits = 1) + ' M '), 
							hjust = just,
							vjust = .5, size = 2) +
		labs(x = NULL, y = TeX(r'(Relative fluorescence $\log_2\frac{F}{F_0}$ (95%CI) )')) +
		scale_x_discrete(labels = parse(text = solute.label(summ$solute))) +
		theme(legend.position = "none", axis.text.x = element_text(angle = 90))
}

dt.aa = dt.bar.summary(aa.450.514, aa.o.water.450.514, aa.conc)
dt.o = dt.bar.summary(rbind(o.450.514, alb.450.514), aa.o.water.450.514, c(o.conc, alb.conc))
dt.salt = dt.bar.summary(salt.450.514, salt.water.450.514, salt.conc)
bar.plot.report(rbind(dt.aa, dt.o, dt.salt))
save.plot(directory('Figures/Report/') + '514 barplots (ex 450).png',
					width = 6.31, height = 2.5, units = 'in')

dt.aa = dt.bar.summary(aa.405.514, aa.o.water.405.514, aa.conc)
dt.o = dt.bar.summary(rbind(o.405.514, alb.405.514), aa.o.water.405.514, c(o.conc, alb.conc))
dt.salt = dt.bar.summary(salt.405.514, salt.water.405.514, salt.conc)
bar.plot.report(rbind(dt.aa, dt.o, dt.salt))
save.plot(directory('Figures/Report/') + '514 barplots (ex 405).png',
					width = 6.31, height = 2.5, units = 'in')


#AWE, Serum, PBS
p1 = emission.plot(ph.pbs.serum.awe[wl.nm <= 540], c('solution', 'ph', 'wl.nm'), 
							aes(wl.nm, color = solution, fill = solution, 
									linetype = signif.as.character(ph, 2)), 
							fill.label = '', color.label = '', linetype = 'pH') +
	annotate('rect', xmin = -Inf, xmax = 514, ymin = 30000, ymax = Inf, fill = 'magenta', alpha = .1) +
	annotate('rect', xmin = -Inf, xmax = 514, ymin = 5000, ymax = 30000, fill = 'yellow', alpha = .1) +
	annotate('rect', xmin = -Inf, xmax = 514, ymin = -Inf, ymax = 5000, fill = 'dark blue', alpha = .1) +
	geom_text(aes(x = 514, y = 30000), label = 'Chronic?', hjust = 1, vjust = 0, size = 2.5, 
						color = 'magenta') +
	geom_text(aes(x = 514, y = 5000), label = 'Healing?', hjust = 1, vjust = 1, size = 2.5, 
					color = 'dark blue')

p2 = emission.plot(ph.pbs.serum.awe[wl.nm <= 540], c('solution', 'ph', 'wl.nm'), 
									aes(wl.nm, color = solution, fill = solution), 
									fill.label = '', color.label = '') +
	facet_wrap(~ fct_rev(paste0('pH ', signif.as.character(ph, 2))), scales = 'free_y', ncol = 1) +
	theme(axis.title.y = element_blank())

p1 + p2 + plot_layout(guides = 'collect', widths = c(5, 3))
save.plot(directory('Figures/Report/') + 'AWE Serum PBS.png',
					width = 6.31, height = 3, units = 'in')

#Resazurin assay optimization----
res.opt.abs.xl = read_excel('Data/R/Resazurin optimization.xlsx', sheet = 'Absorbance')
res.opt.flu.xl = read_excel('Data/R/Resazurin optimization.xlsx', sheet = 'Fluorescence')
res.opt.abs = xl.to.dt(res.opt.abs.xl, 'A')
res.opt.flu = xl.to.dt(res.opt.flu.xl)

res.data = function(dt){
	dt[type == 'Cells', -'type']
}

res.nc = function(dt){
	dt = dt[type == 'Negative control' & dye.M != 0, -'type']
	if('cells.per.ml' %in% colnames(dt)){
		dt = dt[, -'cells.per.ml']
	}
	if('hacc.ppm' %in% colnames(dt)){
		dt = dt[, -'hacc.ppm']
	}
	dt
}

res.pc = function(dt){
	dt = dt[type == 'Positive control', -'type']
	if('cells.per.ml' %in% colnames(dt)){
		dt = dt[, -'cells.per.ml']
	}
	if('hacc.ppm' %in% colnames(dt)){
		dt = dt[, -'hacc.ppm']
	}
	dt
}

#Absorbance method
red.percent.abs = function(dt, h, cells, d){
	data = res.data(dt)[hour == h & cells.per.ml == cells & dye.M == d]
	nc = res.nc(dt)[hour == h & dye.M == d]
	e.o.570 = 80586
	e.o.600 = 117216
	e.r.570 = 155677
	e.r.600 = 14652
	a.570 = data[wl.nm == 570, o.d.]
	a.600 = data[wl.nm == 600, o.d.]
	c.570 = nc[wl.nm == 570, o.d.]
	c.600 = nc[wl.nm == 600, o.d.]
	(e.o.600 * a.570 - e.o.570 * a.600) / (e.r.570 * c.600 - e.r.600 * c.570) * 100
}

#Fluorescence method
red.percent.flu = function(dt, h, d, f){
	nc = res.nc(dt)[hour == h & dye.M == d, rfu]
	pc = res.pc(dt)[hour == h & dye.M == d, rfu]
	(f - nc) / (pc - nc) * 100
}

red.p.dt = function(dt, mode = 'fluorescence'){
	data = res.data(dt)
	result.dt = NULL
	for(h in unique(data$hour)){
		for(cells in unique(data$cells.per.ml)){
			for(d in unique(data$dye.M)){
				if(nrow(data[hour == h & cells.per.ml == cells & dye.M == d]) > 0){
					dt.n = data.table(hour = h, cells.per.ml = cells, dye.M = d,
														red.p = switch(mode,
																					 'absorbance' = red.percent.abs(dt, h, cells, d),
																					 'fluorescence' = red.percent.flu(dt, h, cells, d)$mean))
					result.dt = rbind(result.dt, dt.n)
				}
			}
		}
	}
	remove.na(result.dt)
}
dt = res.opt.flu[wl.nm == 560]
res.opt.red.p.f = 
	res.data(dt)[, .(red.p = mean(red.percent.flu(dt, hour, dye.M, rfu)), mode = 'Fluorescence'),
							keyby = .(hour, cells.per.ml, dye.M)]
dt = res.opt.abs
res.opt.red.p.a = 
	res.data(dt)[wl.nm == 570][, .(red.p = mean(red.percent.abs(dt, hour, cells.per.ml, dye.M)),
																 mode = 'Absorbance'),
							keyby = .(hour, cells.per.ml, dye.M)]
res.opt.red.p = rbind(res.opt.red.p.a, res.opt.red.p.f)

#Faced grid plot
p1 = res.opt.red.p[hour < 8] %>%
	ggplot(aes(hour, red.p, shape = as.factor(cells.per.ml), color = as.factor(cells.per.ml))) +
	geom_point(alpha = .5) + geom_line(alpha = .5) +
	labs(x = 'Incubation time [h]', y = 'Resazurin % reduction', 
			 shape = 'Cells/ml', color = 'Cells/ml') +
	facet_nested('Resazurin concentration [M]' + dye.M ~ mode, scales = 'free') + 
	scale_x_continuous(breaks = 0:7) +
	theme(strip.background.y = element_blank(), strip.text.y = element_blank())

p2 = res.opt.red.p[hour < 8] %>%
	ggplot(aes(cells.per.ml, red.p, color = as.factor(hour))) +
	geom_point(alpha = .5) + geom_line(alpha = .5) +
	labs(x = 'Cells/ml', y = 'Resazurin % reduction', color = 'Time [h]') +
	facet_nested('Resazurin concentration [M]' + dye.M ~ mode, scales = 'free') + 
	scale_x_continuous(breaks = c(10000, 20000, 40000, 80000)) +
	theme(axis.text.x = element_text(angle = 90), axis.text.y = element_blank(), 
				axis.title.y = element_blank(), axis.ticks.y = element_blank())

p1 + p2 + plot_layout(guides = 'collect')
save.plot(directory('Figures/Report/') + 'Resazurin optimization.png',
					width = 6.31, height = 4, units = 'in')

#Regression
dt = res.opt.red.p.f[dye.M == 4e-6]
dt.reg = dt[cells.per.ml %in% c(5000, 20000, 80000)]
dt.reg %>%	ggplot(aes(cells.per.ml, red.p, color = as.factor(hour))) +
	geom_point() +
	geom_smooth(method = 'lm', linewidth = .5, alpha = .5, data = dt.reg) + 
	stat_regline_equation(aes(label =  ..adj.rr.label..),
												label.y = Inf, vjust = 1, data = dt.reg, size = 2.5) +
	stat_cor(aes(label =  paste(..p.label.., sep = "~~~~")),
					 label.y = Inf, vjust = 2.5, data = dt.reg, size = 2.5) +
	labs(x = 'Cells/ml', y = 'Resazurin % reduction') +
	facet_wrap(~ paste0(hour, ' h'), ncol = 3) + 
	scale_x_continuous(breaks = unique(dt.reg$cells.per.ml)) +
	scale_y_continuous(breaks = c(0, 50, 100)) +
	theme(legend.position = 'none', axis.text.x = element_text(angle = 90))
save.plot(directory('Figures/Report/') + 'Resazurin optimization regression.png',
					width = 6.31, height = 3, units = 'in')


#Resazurin assay----
dt.res.xl = read_excel('Data/R/Resazurin assay.xlsx', sheet = 'Fluorescence')
dt.res = xl.to.dt(res.flu.xl)[wl.nm == 560] %>% mutate(inc.hour = hour %% 24, wl.nm = NULL)
dt.res = dt.res[order(z.position, gain, hour, type, hacc.ppm, dye.M)]

#Data filtered for gain and z-position
dt = dt.res[gain == 70 & z.position == 15000] %>% dplyr::select(-gain, -z.position)

#Data filtered for wells without treatment
dt = dt[hour < 40 | column != 4 | row %in% c('A', 'B', 'C')]

#Data filtered for incubation time
dt = dt[inc.hour == 4] %>% dplyr::select(-inc.hour)

#Data filtered for missing dye
dt.m = dt[dye.M == 0 | rfu >= 1000]

#Data without outliers
dt.ol = dt.m %>% remove.outliters(c('hour', 'type', 'hacc.ppm', 'dye.M'), 'rfu')

#Diagnostic plots
#Gain comparison
gain.comp = dt.res[hour %in% c(4, 25, 28, 49, 52)]
gain.comp$type = factor(gain.comp$type, 
												levels = c('Cells', 'Negative control', 'Positive control'),
												labels = c('Cells', 'N. control', 'P. control'))
gain.comp = cbind(gain.comp[gain == 100] %>% mutate(g.100 = rfu) %>% dplyr::select(-rfu), 
									(gain.comp[gain == 70] %>% mutate(g.70 = rfu)) %>% dplyr::select(g.70))
gain.comp.low = gain.comp[g.70 < 50000]
gain.comp.high = gain.comp[g.70 >= 50000]
p1 = gain.comp %>%
	ggplot(aes(g.70, g.100, color = type)) +
	geom_point(alpha = .2) +
	geom_smooth(method = 'lm', linewidth = .3, color = 'purple') +
	stat_cor(aes(label =  paste(..rr.label.., ..p.label.., sep = "~`, `~")),
					 color = 'purple', data = gain.comp, label.x = 40000, label.y = 1250000, size = 2) +
	
	geom_smooth(method = 'lm', linewidth = .3, color = 'red', data = gain.comp.low) +
	stat_cor(aes(label =  paste(..rr.label.., ..p.label.., sep = "~`, `~")),
					 color = 'red', data = gain.comp.low, label.y = 700000, size = 2) +
	
	geom_smooth(method = 'lm', linewidth = .3, color = 'blue', data = gain.comp.high) +
	stat_cor(aes(label =  paste(..rr.label.., ..p.label.., sep = "~`, `~")), color = 'blue', 
					 data = gain.comp.high, label.x = 90000, label.y = 1000000, size = 2) +
	
	labs(x = 'Gain 70 fluorescence [RFU]', y = 'Gain 100 fluorescence [RFU]', color = '') +
	scale_x_continuous(labels = function(x) format(x, scientific = TRUE)) +
	scale_y_continuous(labels = function(x) format(x, scientific = TRUE))

#Z-position comparison
z.comp = res.flu[hour %in% c(25, 28, 49, 52)]
z.comp = cbind(z.comp[z.position == 15000] %>% mutate(z.manual = rfu) %>%	dplyr::select(-rfu), 
							 z.comp[z.position == 'auto'] %>% mutate(z.auto = rfu) %>% dplyr::select(z.auto))
z.comp = z.comp[order(z.auto)]
z.comp.low = z.comp[z.auto < 750000]
z.comp.high = z.comp[z.auto >= 750000]
p2 = z.comp %>%
	ggplot(aes(z.auto, z.manual, color = as.factor(gain))) +
	geom_point(alpha = .1) +
	geom_smooth(method = 'lm', linewidth = .3, color = 'purple') +
	stat_cor(aes(label =  paste(..rr.label.., ..p.label.., sep = "~`, `~")),
					 color = 'purple', label.x = 500000, label.y = 1100000, size = 2) +
	geom_smooth(method = 'lm', linewidth = .3, color = 'red', data = z.comp.low) +
	stat_cor(aes(label =  paste(..rr.label.., ..p.label.., sep = "~`, `~")),
					 color = 'red', data = z.comp.low, label.y = 600000, size = 2) +
	geom_smooth(method = 'lm', linewidth = .3, color = 'blue', data = z.comp.high) +
	stat_cor(aes(label =  paste(..rr.label.., ..p.label.., sep = "~`, `~")), color = 'blue',
					 data = z.comp.high, label.x = 1000000, label.y = 960000, size = 2) +
	labs(x = 'Auto-Z fluorescence [RFU]', y = 'Manual-Z fluorescence [RFU]', 
			 color = 'Gain') +
	scale_x_continuous(labels = function(x) format(x, scientific = TRUE)) +
	scale_y_continuous(labels = function(x) format(x, scientific = TRUE))

p1 + p2 + plot_annotation(tag_levels = 'A') & theme(legend.position = 'bottom')
save.plot(directory('Figures/Report/') + 'Resazurin gain z.png',
					width = 6.31, height = 3.7, units = 'in')

#Combined data
dt.all = rbind(dt %>% mutate(dat = 'None'),
							 dt.m %>% mutate(dat = 'No dye'),
							 dt.ol %>% mutate(dat = 'All'))
dt.all$dat = factor(dt.all$dat, levels = c('None', 'No dye', 'All'))
#Raw plot
ctrl.lbl = function(t) str_replace(
	str_replace(t, 'Positive control', 'PC'), 'Negative control', 'NC')
dye.lbl = function(d) str_replace(str_replace(round(d * 1e6), '4', 'Dye'), '0', '-')

raw.plot = function(dt){
	dt %>% ggplot(aes(hour, rfu)) +
		facet_nested('Removed outliers' + dat ~ ctrl.lbl(type) + round(hacc.ppm) + dye.lbl(dye.M)) +
		geom_point(alpha = .1) +
		geom_boxplot(aes(group = cut_width(hour, 1)), outlier.color = 'red', alpha = .1) +
		labs(x = 'Time since treatment [h]', y = 'Fluorescence [RFU]') +
		scale_x_continuous(breaks = unique(dt$hour)) +
		scale_y_continuous(labels = function(x) format(x, scientific = TRUE, digit = 1))
}
raw.plot(dt.all)
save.plot(directory('Figures/Report/') + 'Resazurin raw plot.png',
					width = 6.31, height = 3.7, units = 'in')

#Missing dye (RFU < 1000)
dt.nc = dt[type == 'Negative control' | type == 'Cells']
dt.nc$f = factor(dt.nc$type, levels = c('Negative control', 'Cells'))
dt.nc$d = sapply(dt.nc$dye.M, function(x) ifelse(x == 0, 'no dye', 'dye'))
dt.nc$d = factor(dt.nc$d, levels = c('no dye', 'dye'))
nc.99.9.high = ci(dt.nc[dye.M == 0]$rfu, .999)$high
p.ol = dt.nc %>%
	ggplot(aes(as.factor(hour), rfu)) +
	geom_boxplot(outlier.color = 'red', outlier.alpha = .2) +
	facet_nested(. ~ f + d) +
	geom_abline(aes(slope = 0, intercept = 1000), linetype = 'dotted') +
	geom_text(x = -Inf, y = 1000, hjust = 0, vjust = 0, size = 2.5, label = '1000') +
	labs(x = 'Time since treatment [h]', y = 'Fluorescence [RFU]')

#Control trends
control.trend.plot = function(dt){
	dt = dt[type %in% c('Negative control', 'Positive control')]
	summ = dt[, .(rfu = mean(rfu), hw = ci(rfu)$hw, min = min(rfu)), by = .(dye.M, type, hour, dat)]
	dt %>%
		ggplot(aes(as.factor(hour), rfu)) +
		geom_boxplot(outlier.color = 'red', outlier.alpha = .2, varwidth = T) +
		geom_errorbar(data = summ, aes(ymin = rfu - hw, ymax = rfu + hw), width = .4, alpha = .2) +
		geom_text(data = summ, aes(y = min, label = round(rfu) + ' ± ' + round(hw)),  
							vjust = 2, hjust =.5,	size = 1.5, alpha = .5) +
		facet_nested('Removed outliers' + dat ~ type + dye.lbl(dye.M)) +
		geom_pwc(aes(group = hour), method = 't_test', tip.length = 0, label.size = 2.5,
						 p.adjust.method = 'BH', bracket.nudge.y = .2) +
		stat_kruskal_test(size = 3, vjust = -1, hjust = .5, label.x = 2, label.y = dt[, max(rfu)]) +
		labs(x = 'Time since treatment [h]', y = 'Fluorescence [RFU]')
}
control.trend.plot(dt.all)
save.plot(directory('Figures/Report/') + 'Resazurin control trend.png',
					width = 6.31, height = 6.7, units = 'in')

#Cell trends
cell.trend.plot = function(dt){
	dt = dt[type == 'Cells']
	dt %>%
		ggplot(aes(as.factor(hour), rfu)) +
		geom_boxplot(outlier.color = 'red', outlier.alpha = .2, varwidth = T) +
		geom_smooth(method = 'lm', formula = rfu ~ exp(hour)) +
		facet_nested('Removed outliers' + dat ~ 'Nanoparticle concentration [ppm]' + round(hacc.ppm)) +
		geom_pwc(aes(group = hour), method = 't_test', tip.length = 0, label.size = 2,
						 p.adjust.method = 'BH', bracket.nudge.y = .2) +
		stat_kruskal_test(label = 'K-W\n{p.format}',
										size = 2, vjust = 0, hjust = .5, label.x = 2, label.y = dt[, max(rfu)]) +
		labs(x = 'Time since treatment [h]', y = 'Fluorescence [RFU]')
}
cell.trend.plot(dt.all)
save.plot(directory('Figures/Report/') + 'Resazurin treatment trend.png',
					width = 6.31, height = 6, units = 'in')


#Fluorescence method
red.percent = function(dt, h, hacc){
	f = res.data(dt)[hour == h & hacc.ppm == hacc]$rfu
	nc = res.nc(dt)[hour == h]$rfu
	pc = res.pc(dt)[hour <= 28]$rfu
	error.propagation('(f - nc) / (pc - nc) * 100', f, nc, pc)
}

red.p.summary = function(dt){
	res.data(dt)[, red.percent(dt, hour, hacc.ppm),	keyby = .(hour, hacc.ppm)]
}

summ = rbind(red.p.summary(dt) %>% mutate(dat = 'None'), 
						 red.p.summary(dt.m) %>% mutate(dat = 'No dye'), 
						 red.p.summary(dt.ol) %>% mutate(dat = 'All'))
summ$dat = factor(summ$dat, levels = c('None', 'No dye', 'All'))

res.colors = rainbow(length(unique(summ$hacc.ppm)))
names(res.colors) = as.character(round(unique(summ$hacc.ppm)))

#Faced wrap plot
#by ppm
red.p.facet.ppm = function(dt){
	dt %>%
		ggplot(aes(hour, mean, color = as.factor(round(hacc.ppm)))) +
		geom_point(alpha = .5, size = 1) + geom_line(alpha = .5) +
		geom_ribbon(aes(ymin = low, ymax = high, fill = as.factor(round(hacc.ppm))), 
								alpha = 0.2, color = NA) +
		geom_errorbar(aes(ymin = low, ymax = high), width = 0, alpha = 0.5) +
		facet_nested('Removed outliers' + dat ~ 'Treatment concentration [ppm]' + round(hacc.ppm)) +
		labs(x = 'Time since treatment [h]', y = 'Resazurin % reduction (95%CI)') +
		scale_x_continuous(breaks = unique(dt$hour)) +
		scale_color_manual(values = res.colors) +
		scale_fill_manual(values = res.colors) +
		theme(legend.position = 'none', axis.text.x = element_text(angle = 90),
					strip.background.y = element_blank(), strip.text.y = element_blank())
}

#by hour
red.p.facet.hour = function(dt){
	ctrl = dt[hacc.ppm == 0]
	dt = dt[hacc.ppm != 0]
	ppm = unique(dt$hacc.ppm)
	dt %>%
		ggplot(aes(hacc.ppm, mean, color = hour, fill = hour)) +
		geom_point(alpha = .5, size = 1) + geom_line(alpha = .5) +
		geom_ribbon(aes(ymin = low, ymax = high), alpha = 0.2, color = NA) +
		geom_errorbar(aes(ymin = low, ymax = high), width = 0, alpha = 0.5) +
		geom_hline(aes(yintercept = mean), alpha = .5, data = ctrl) +
		geom_rect(aes(xmin = 0, xmax = Inf, ymin = low, ymax = high), 
							alpha = .1, color = NA, fill = 'black', data = ctrl) +
		facet_nested('Removed outliers' + dat ~ 'Time since treatment [h]' + hour) +
		labs(x = 'Treatment concentration [ppm]', y = 'Resazurin % reduction (95%CI)') +
		scale_x_continuous(breaks = unique(round(ppm)), trans = log_trans(base = 2)) +
		theme(axis.text.x = element_text(angle = 90), legend.position = 'none',
					axis.title.y = element_blank(), axis.text.y = element_blank(), 
					axis.ticks.y = element_blank())
}

red.p.facet.ppm(summ) +
red.p.facet.hour(summ) & 
	theme(strip.text = element_text(size = 6), axis.text.x = element_text(size = 6))
save.plot(directory('Figures/Report/') + 'Resazurin by ppm and hour.png',
					width = 6.31, height = 3, units = 'in')
