# fig 5 c

windows()
# tiff("Fig5d.tiff", height = 20, width = 34, units = 'cm',compression = "lzw", res = 300)

data = c(.05,.05,0.008,.1677,.001)
model = c(.048606,.048296,0.0083882,0.172991,0.001005)


par(mar=c(7,5,5,3))
barplot(data,
	main = 'Comparison between model and data', cex.main = 3,
	ylim=c(0,.2),cex.axis = 2,
	)
strg = expression(frac('PI(4,5)P'[2],'PItotal'))
x_labels = c(
		expression(frac('PI(4,5)P'[2],'PItotal')),
		expression(frac('PI(4)P','PItotal')),
		expression(frac('DAG','PItotal')),
		expression(frac('PAtotal','PItotal')),
		expression(frac('CDPDAG','PItotal'))
		)
text(x = c(.7,1.9,3.1,4.3,5.5),y = rep(-0.01,5), labels = x_labels, srt = 0, pos = 1, xpd = TRUE,cex=2)

segments(x0 = c(.7,1.9,3.1,4.3,5.5)-.3, y0 = model, x1 = c(.7,1.9,3.1,4.3,5.5)+.3, y1 = model,col = 'blue', lwd = 5)

legend(
	x=.2, y=.21,
	legend=c('Data'),
	text.col=c('black'),
	pch=c(15),
	col=c('grey'),
	bty = "n",
	cex=3
	)
legend(
	x=.1, y=.19,
	legend=c('Model'),
	text.col=c('black'),
	lty=1,lwd=5,seg.len=.5,
	col=c('blue'),
	bty = "n",
	cex=3,pt.cex=3
	)


# dev.off()