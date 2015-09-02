#STARTTIME
2009		iYear
07		iMonth
22		iDay
00		iHour
00		iMinute
00		iSecond

#STOP
100.0		TimeMax


#NGDC_INDICES
Indices.dat	NameNgdcFile

#MHD_INDICES
imf.dat		UpstreamFiles

#SMOOTH
F		UseSmooth
300.0		SmoothWindow

#BMODEL
t04		NameModel

#IEMODEL
T		UseWeimer

#SAVEPLOT
10.0		DtSavePlot
T		DoSaveFlux
F		DoSaveDrifts
F		UseSeparatePlotFiles

#SAVELOG
10.0		DtLogOut

#INITIALF2
F		IsEmptyInitial
F		IsGmInitial
T		IsDataInitial

#TYPEBOUNDARY
Ellipse		TypeBoundary

#PLASMASHEET
T		UseYoungEtAl
T		UseBoundaryEbihara

#RESTART
F		IsRestart

#SAVERESTART
100.0		DtSaveRestart

#TIMESIMULATION
0.0		Time

