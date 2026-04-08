#
#	                  Makefile
#	                  fuer Programm mg40
#	                  erzeugt mit Makefile-Generator
#
#  dated:  Mi Jul 20 13:41:03 CEST 2005
#
.SUFFIXES:  .f .o .c .src .exe .out .z67 .i .l .list
#
#	                  Welche Shell fuer Kommandos  
SHELL = /bin/sh
#
#	                  Namen des Programms definieren
NAME= ../mg40
#
#include mglet.altix
include mglet.pc
#######################################################
#
#	                  Datei auf die die Ladeliste
#	                  ausgegeben wird
MAPFILE= mg40.map,part
#
OBJECTS= \
addmpi.o  \
addregion.o  \
adjdpdx.o \
aeqgrd.o  \
amult.o  \
amult2.o  \
apfil.o  \
aream.o  \
asdche.o  \
asdfun.o  \
asdval.o  \
baukor.o  \
bbacmg.o  \
bbacon.o  \
bbafix.o \
bbakonv.o  \
bbanos.o  \
bbaop1.o  \
bbaop2.o  \
bbaop3.o  \
bbapar.o  \
bbaper.o  \
bbasca.o  \
bboblo.o  \
bbofix.o  \
bbonos.o  \
bboop1.o  \
bbopar.o  \
bbosca.o  \
bbosli.o  \
bbotmg.o  \
bcubmg.o  \
bcunos.o  \
bcusca.o  \
bcusli.o  \
bf1cub.o  \
bfrblo.o  \
bbablo.o  \
briblo.o  \
bleblo.o  \
bfrcon.o  \
bfrfix.o  \
bfrnos.o  \
bfrop1.o  \
bfromg.o  \
bfrpar.o  \
bfrper.o  \
bfrsca.o  \
berrhs.o  \
bl_diagnose.o  \
blecon.o  \
blenos.o  \
bleop1.o  \
bleop2.o  \
blepar.o  \
bleper.o  \
blesca.o  \
blesli.o  \
blftmg.o  \
bodycheck.o  \
bodyst.o  \
boflu.o  \
borec.o  \
boundmg.o  \
boundscainter.o  \
bparmg.o  \
bpart.o  \
bpart_body.o  \
bphi0.o  \
brenner.o \
brgtmg.o  \
bricon.o  \
brifix.o  \
brinos.o  \
briop1.o  \
bripar.o  \
briper.o  \
brisca.o  \
brisli.o  \
btonos.o  \
btoop1.o  \
btoop2.o  \
btoop3.o  \
btopar.o  \
btopmg.o  \
btosca.o  \
btosli.o  \
cal_rr.o  \
cal_rrrr.o  \
coefdiv.o  \
cotofine3d.o  \
calcoefd.o  \
coefkon.o  \
coefscal.o  \
cop3dzero.o  \
calcoefdx.o  \
calcoeff.o  \
calcoeffx.o  \
caldpdx.o  \
calpsfak.o  \
cap1252.o  \
cap1262.o  \
channeld.o  \
coeffop3.o  \
coeffop32.o  \
compose3d.o  \
composefield.o  \
composegrd.o  \
conbacmg.o  \
conbapar.o  \
conbotmg.o  \
confromg.o  \
conlftmg.o  \
connectmg.o  \
confrpar.o  \
conrgtmg.o  \
contopar.o  \
contopmg.o  \
cop3d.o  \
copfib.o  \
copfic.o  \
copyregion.o  \
cousin.o  \
cubinf.o  \
d20.o  \
deconvolv2.o  \
deconvolvsca.o  \
deibi.o  \
deici.o  \
deobi.o  \
deoci.o  \
dfdx.o  \
dib.o  \
dibca.o  \
dibca1.o  \
dibhead.o  \
dic.o  \
dichead.o  \
dicrec.o  \
diffcoefkobs0.o    \
digrid.o  \
dimlo1.o  \
dimlos.o  \
dissipg.o  \
dissipga.o  \
divcal.o  \
dmixle.o  \
dob.o  \
dobca.o  \
dobca1.o  \
dobhead.o  \
dobrec.o  \
doc.o  \
dochead.o  \
doenst.o  \
dogrid.o  \
dphi0.o  \
dtmxcal.o  \
dumpfi.o  \
efvisc.o  \
enerfg.o  \
enerfs.o  \
enstro.o  \
ernorm.o  \
errr.o  \
exchan.o  \
exchansca.o  \
expan1.o  \
extendk.o \
fderfouvz.o  \
fderfouwy.o  \
fderfouwyper.o  \
fderfoux.o  \
fderfouxper.o  \
fderfosca.o  \
fderfovwx.o  \
fderfovwxper.o  \
fderfovy.o  \
fderfovyper.o  \
fderfowz.o  \
fftpack.o  \
findposition.o  \
filter_explizit.o  \
flalfa.o  \
flgama.o  \
fluctuations.o  \
fmgout.o  \
fsensi.o  \
fsensj.o  \
fsensk.o  \
funduc.o  \
galirt.o  \
getderivatives.o  \
getvelocities.o  \
getscalar.o  \
grdbko.o  \
grdche.o  \
grdchk.o  \
grdctof.o  \
grdfit.o  \
grdfmi.o  \
grdftoc.o  \
grdkon.o  \
grdkor.o  \
grdout.o  \
grdrec.o  \
gridae.o  \
gsit.o  \
giteig.o  \
helici.o  \
hrelom.o  \
htmles.o  \
ibfield.o  \
icobody.o  \
icobound.o  \
icogrdcon.o  \
icogrddef.o  \
icogrdpro.o  \
icolevel.o  \
icomgrid.o  \
icomgvp.o  \
icophyspar.o  \
inigrid.o  \
inislice.o  \
initsca.o  \
initorrsommer.o  \
init_part.o  \
init_filter_periodic.o  \
intercoef1.o  \
intercoef2.o  \
intercoef3.o  \
intercoef4.o  \
intercoef5.o  \
intercoefkobs0.o\
intercoefku3.o \
interpolate1.o  \
interpolate2.o  \
interpolate3.o  \
interpolatedj.o  \
interpolatedk.o  \
interpolatekj.o  \
interpolatekk.o  \
interpolatep.o  \
interpolateuvz.o  \
interpolateuwy.o  \
interpolatesca.o  \
interpolateux.o  \
interpolateux3.o \
interpolatevwx.o  \
interpolatevy.o  \
interpolatewz.o  \
interpolphi.o  \
interuwyper.o  \
interuxper.o  \
intervwxper.o  \
intervyper.o  \
interpolate_3d_cell2.o \
interpolate_3d_cell3.o \
interpolate_3d_cell6.o \
intlin.o  \
intst.o  \
intsta.o  \
itinf.o  \
itsample.o  \
kjidco.o  \
kjieco.o  \
ko2hom.o  \
ko2pkt.o  \
ko2var.o  \
kstprt.o  \
lesconst.o  \
linctl.o  \
lininf.o  \
linout.o  \
listi6.o  \
mamili.o  \
mgbasb.o  \
mgbasbsca.o  \
mgbftc.o  \
mgbloindset.o  \
mgboflu.o  \
mgconinf.o  \
mgcover.o  \
mgctof.o  \
mgdima.o  \
mgdims.o  \
mgdpb.o  \
mgftoc.o  \
mggrdgen.o  \
mgnbrbuf.o  \
mgnbrcheck.o  \
mgnbrset.o  \
mgoverlap.o  \
mgparcheck.o  \
mgparset.o  \
mgpcorr.o  \
mgpoina.o  \
mgpoint.o  \
mgpoisit.o \
mgpoisc1.o  \
mgpoisl1.o  \
mgpsdir.o  \
mgsvflu.o  \
mgtstpart.o \
mgvpc1.o  \
mgvpc2.o  \
mgvpc3.o  \
mgvpc4.o  \
mgvpc5.o  \
mgvpit.o  \
mgvpl1.o  \
mgvpset.o  \
mlet.o  \
msgvaropt.o \
nextgrid.o  \
omega.o  \
optfreq.o \
output.o  \
outres.o  \
orrsommer.o  \
parbacmg.o  \
parbotmg.o  \
parfromg.o  \
parlftmg.o  \
parrgtmg.o  \
parset.o  \
partdiss.o  \
particle_stress.o \
partopmg.o  \
passpart.o  \
perturb.o  \
phiadd.o  \
phididj.o  \
phifl2.o  \
phifla.o  \
phiflu.o  \
phimlt2.o \
phimlt3.o \
phimlt4.o \
phirmi.o  \
phiske.o  \
plevel.o  \
pnivea.o  \
posrec.o  \
pr1lin.o  \
pr3b.o  \
preproc_field.o \
preproc_field_scalar.o \
preproc_field_wcomp.o \
preproc_pressure.o \
pr3v.o  \
printe.o  \
printvel.o  \
prle3e.o  \
prle3l.o  \
prle3m.o  \
prolong1.o  \
prolong2.o  \
random.o  \
random_functions.o \
random_numbers.o \
randwert.o  \
randwertan.o  \
randwertand.o  \
randwertd.o  \
randwertdw.o  \
randwertij.o  \
randwertijan.o  \
randwertijand.o  \
randwertijandw.o  \
randwertijd.o  \
randwertijdw.o  \
randwertijw.o  \
randwertw.o  \
read_particles.o \
read_particles_ini.o  \
readsliced.o  \
recout.o  \
rezipd.o  \
rkjiae.o  \
rmcomment.o  \
sbbc12.o  \
sbcoun.o  \
sbcub.o  \
sbhemi.o  \
sbzyl.o  \
searchindex.o \
sel0.o  \
sel1.o  \
sel2.o  \
sel3.o  \
sel4.o  \
sel5.o  \
sel6.o  \
sel7.o  \
sel8.o  \
sel11.o \
sel12.o \
sel13.o \
sel14.o \
sel15.o \
sel17.o \
sel20.o \
selaco.o  \
selauf.o  \
seleca.o  \
seleci.o  \
setbackold.o  \
setbloindex.o  \
setcobone.o  \
setcobound.o  \
setcolevel.o  \
setcomgrid.o  \
setgeovp.o  \
setgrd.o  \
setgrdpro.o  \
setid8.o  \
setidx.o  \
setkon.o  \
setmpi.o  \
setrandom.o  \
setref.o  \
sets.o  \
setslice.o  \
setst1.o  \
setsta.o  \
setthreads.o  \
sipiter.o  \
siplu.o  \
skonle.o  \
slice1d.o  \
slice3d.o  \
slicefield.o  \
slicegeo.o  \
slicegrd.o  \
smooth_stress.o \
srmsgs.o  \
stabed.o  \
startcontrol.o  \
stmim1.o  \
stmimp.o  \
stmnew.o  \
stmnli.o  \
stmnrm.o  \
strcheck.o  \
stress.o \
strles.o  \
strset.o  \
sum2ar.o  \
sumaan.o  \
sumst1.o  \
sumsta.o  \
sv3d.o  \
sveipr.o  \
sveiprsca.o  \
svflu.o  \
svle1.o  \
svle1sca.o  \
svrec.o  \
swcle3d.o  \
swclesca.o  \
taufg.o  \
taufm.o  \
taufs.o  \
taunn.o  \
thomasi.o  \
thomasj.o  \
thomask.o  \
trizyk.o  \
trizyk3d.o  \
trizyk3dj.o  \
tst1g.o  \
tst3rk.o  \
tstle2.o  \
tstle4.o  \
tstorient.o  \
tstpos.o  \
tstsca.o  \
tstsca4.o  \
ubulk.o  \
vplesc.o  \
vprbc.o  \
vshift.o  \
wallkobdiffs0a.o  \
wallkobdiffs0e.o   \
wifak.o  \
wifaksca.o  \
wnsint.o  \
wrigeo.o  \
wrirec.o  \
write_dupart.o \
write_particles.o \
write_particles_fin.o  \
write_passpart.o  \
write3d.o  \
writesliced.o  \
wrivaropt.o \
wssint.o  \
xtractxy.o  \
xtractyz.o  \
xtractyzf.o  \
zbrent.o   \
dyn.o \
inc.o \
inch.o \
leon.o \
model.o \
tauinc.o \
conrgtmgbp.o \
confromgbp.o \
conbacmgbp.o \
conlftmgbp.o \
briconbp.o \
bleconbp.o \
bfrconbp.o \
bbaconbp.o \
inigridgeo.o \
inigridvar.o \
blockbp.o  \
filter_periodicx.o  \
filter_periodicy.o \
printdata.o 



#FFT_OBJECTS = \
radb2.o \
radb3.o \
radb4.o \
radb5.o \
radbg.o \
radf2.o \
radf3.o \
radf4.o \
radf5.o \
radfg.o \
rfftb.o \
rfftb1.o \
rfftf.o \
rfftf1.o \
rffti.o \
rffti1.o
# filter_periodicx.o  \
# filter_periodicy.o  \

MAIN=\
mlet.o

MGLET_PRE=\
mglet_pre.o

#	                  Compilier-Aufruf
#
mg40.exe:  $(OBJECTS) $(FFT_OBJECTS) mg40.def ./qpack/libQPack.a
	$(F77) $(OBJECTS) $(FFT_OBJECTS) $(LOPTIONS) -o $(NAME).exe $(LIBRARY)
#
addmpi.o: addmpi.src mg40.def cogrdpro.h colevel.h compi.h mgpar.h
	cat mg40.def addmpi.src |$(CPP) $(CPPFLAGS)  > addmpi.f
	 $(F77) $(OPTIONS) $(FFLAGS)  addmpi.f  
addregion.o: addregion.src mg40.def
	cat mg40.def addregion.src |$(CPP) $(CPPFLAGS)  > addregion.f
	 $(F77) $(OPTIONS) $(FFLAGS)  addregion.f  
adjdpdx.o: adjdpdx.src mg40.def
	cat mg40.def adjdpdx.src |$(CPP) $(CPPFLAGS)  > adjdpdx.f
	 $(F77) $(OPTIONS) $(FFLAGS)  adjdpdx.f  
aeqgrd.o: aeqgrd.src mg40.def colevel.h compi.h konsta.h mgpar.h
	cat mg40.def aeqgrd.src |$(CPP) $(CPPFLAGS)  > aeqgrd.f
	 $(F77) $(OPTIONS) $(FFLAGS)  aeqgrd.f  
amult.o: amult.src mg40.def
	cat mg40.def amult.src |$(CPP) $(CPPFLAGS)  > amult.f
	 $(F77) $(OPTIONS) $(FFLAGS)  amult.f  
amult2.o: amult2.src mg40.def
	cat mg40.def amult2.src |$(CPP) $(CPPFLAGS)  > amult2.f
	 $(F77) $(OPTIONS) $(FFLAGS)  amult2.f  
apfil.o: apfil.src mg40.def
	cat mg40.def apfil.src |$(CPP) $(CPPFLAGS)  > apfil.f
	 $(F77) $(OPTIONS) $(FFLAGS)  apfil.f
aream.o: aream.src mg40.def
	cat mg40.def aream.src |$(CPP) $(CPPFLAGS)  > aream.f
	 $(F77) $(OPTIONS) $(FFLAGS)  aream.f  
asdche.o: asdche.src mg40.def
	cat mg40.def asdche.src |$(CPP) $(CPPFLAGS)  > asdche.f
	 $(F77) $(OPTIONS) $(FFLAGS)  asdche.f  
asdfun.o: asdfun.src mg40.def konsta.h
	cat mg40.def asdfun.src |$(CPP) $(CPPFLAGS)  > asdfun.f
	 $(F77) $(OPTIONS) $(FFLAGS)  asdfun.f  
asdval.o: asdval.src mg40.def
	cat mg40.def asdval.src |$(CPP) $(CPPFLAGS)  > asdval.f
	 $(F77) $(OPTIONS) $(FFLAGS)  asdval.f  
baukor.o: baukor.src mg40.def
	cat mg40.def baukor.src |$(CPP) $(CPPFLAGS)  > baukor.f
	 $(F77) $(OPTIONS) $(FFLAGS)  baukor.f  
bbacmg.o: bbacmg.src mg40.def cobound.h comgrid.h mgpar.h
	cat mg40.def bbacmg.src |$(CPP) $(CPPFLAGS)  > bbacmg.f
	 $(F77) $(OPTIONS) $(FFLAGS)  bbacmg.f  
bbacon.o: bbacon.src mg40.def colevel.h compi.h mgpar.h
	cat mg40.def bbacon.src |$(CPP) $(CPPFLAGS)  > bbacon.f
	 $(F77) $(OPTIONS) $(FFLAGS)  bbacon.f  
bbafix.o: bbafix.src mg40.def colevel.h compi.h mgpar.h
	cat mg40.def bbafix.src |$(CPP) $(CPPFLAGS)  > bbafix.f
	$(F77) $(OPTIONS) $(FFLAGS)  bbafix.f
bbakonv.o: bbakonv.src mg40.def costrles.h
	cat mg40.def bbakonv.src |$(CPP) $(CPPFLAGS)  > bbakonv.f
	$(F77) $(OPTIONS) $(FFLAGS)  bbakonv.f  
bbanos.o: bbanos.src mg40.def konsta.h
	cat mg40.def bbanos.src |$(CPP) $(CPPFLAGS)  > bbanos.f
	$(F77) $(OPTIONS) $(FFLAGS)  bbanos.f  
bbaop1.o: bbaop1.src mg40.def cophyspar.h
	cat mg40.def bbaop1.src |$(CPP) $(CPPFLAGS)  > bbaop1.f
	 $(F77) $(OPTIONS) $(FFLAGS)  bbaop1.f  
bbaop2.o: bbaop2.src mg40.def
	cat mg40.def bbaop2.src |$(CPP) $(CPPFLAGS)  > bbaop2.f
	 $(F77) $(OPTIONS) $(FFLAGS)  bbaop2.f  
bbaop3.o: bbaop3.src mg40.def
	cat mg40.def bbaop3.src |$(CPP) $(CPPFLAGS)  > bbaop3.f
	 $(F77) $(OPTIONS) $(FFLAGS)  bbaop3.f  
bbapar.o: bbapar.src mg40.def
	cat mg40.def bbapar.src |$(CPP) $(CPPFLAGS)  > bbapar.f
	 $(F77) $(OPTIONS) $(FFLAGS)  bbapar.f  
bbaper.o: bbaper.src mg40.def
	cat mg40.def bbaper.src |$(CPP) $(CPPFLAGS)  > bbaper.f
	 $(F77) $(OPTIONS) $(FFLAGS)  bbaper.f  
bbasca.o: bbasca.src mg40.def konsta.h
	cat mg40.def bbasca.src |$(CPP) $(CPPFLAGS)  > bbasca.f
	$(F77) $(OPTIONS) $(FFLAGS)  bbasca.f  	 
bboblo.o: bboblo.src mg40.def cophyspar.h konsta.h
	cat mg40.def bboblo.src |$(CPP) $(CPPFLAGS)  > bboblo.f
	 $(F77) $(OPTIONS) $(FFLAGS)  bboblo.f  
bbofix.o: bbofix.src mg40.def
	cat mg40.def bbofix.src |$(CPP) $(CPPFLAGS)  > bbofix.f
	 $(F77) $(OPTIONS) $(FFLAGS)  bbofix.f  
bbonos.o: bbonos.src mg40.def konsta.h
	cat mg40.def bbonos.src |$(CPP) $(CPPFLAGS)  > bbonos.f
	 $(F77) $(OPTIONS) $(FFLAGS)  bbonos.f  
bboop1.o: bboop1.src mg40.def colevel.h compi.h mgpar.h
	cat mg40.def bboop1.src |$(CPP) $(CPPFLAGS)  > bboop1.f
	 $(F77) $(OPTIONS) $(FFLAGS)  bboop1.f  
bbopar.o: bbopar.src mg40.def colevel.h compi.h mgpar.h
	cat mg40.def bbopar.src |$(CPP) $(CPPFLAGS)  > bbopar.f
	 $(F77) $(OPTIONS) $(FFLAGS)  bbopar.f 
bbosca.o: bbosca.src mg40.def konsta.h
	cat mg40.def bbosca.src |$(CPP) $(CPPFLAGS)  > bbosca.f
	 $(F77) $(OPTIONS) $(FFLAGS)  bbosca.f  	  
bbosli.o: bbosli.src mg40.def konsta.h
	cat mg40.def bbosli.src |$(CPP) $(CPPFLAGS)  > bbosli.f
	 $(F77) $(OPTIONS) $(FFLAGS)  bbosli.f  
bbotmg.o: bbotmg.src mg40.def cobound.h comgrid.h mgpar.h
	cat mg40.def bbotmg.src |$(CPP) $(CPPFLAGS)  > bbotmg.f
	 $(F77) $(OPTIONS) $(FFLAGS)  bbotmg.f  
bcubmg.o: bcubmg.src mg40.def cobound.h comgrid.h mgpar.h
	cat mg40.def bcubmg.src |$(CPP) $(CPPFLAGS)  > bcubmg.f
	 $(F77) $(OPTIONS) $(FFLAGS)  bcubmg.f  
bcunos.o: bcunos.src mg40.def konsta.h
	cat mg40.def bcunos.src |$(CPP) $(CPPFLAGS)  > bcunos.f
	 $(F77) $(OPTIONS) $(FFLAGS)  bcunos.f  
bcusca.o: bcusca.src mg40.def konsta.h
	cat mg40.def bcusca.src |$(CPP) $(CPPFLAGS)  > bcusca.f
	 $(F77) $(OPTIONS) $(FFLAGS)  bcusca.f  	 
bcusli.o: bcusli.src mg40.def konsta.h
	cat mg40.def bcusli.src |$(CPP) $(CPPFLAGS)  > bcusli.f
	 $(F77) $(OPTIONS) $(FFLAGS)  bcusli.f  
bf1cub.o: bf1cub.src mg40.def
	cat mg40.def bf1cub.src |$(CPP) $(CPPFLAGS)  > bf1cub.f
	 $(F77) $(OPTIONS) $(FFLAGS)  bf1cub.f  
bfrblo.o: bfrblo.src mg40.def cophyspar.h
	cat mg40.def bfrblo.src |$(CPP) $(CPPFLAGS)  > bfrblo.f
	 $(F77) $(OPTIONS) $(FFLAGS)  bfrblo.f  
bbablo.o: bbablo.src mg40.def cophyspar.h
	cat mg40.def bbablo.src |$(CPP) $(CPPFLAGS)  > bbablo.f
	 $(F77) $(OPTIONS) $(FFLAGS)  bbablo.f
briblo.o: briblo.src mg40.def cophyspar.h
	cat mg40.def briblo.src |$(CPP) $(CPPFLAGS)  > briblo.f
	 $(F77) $(OPTIONS) $(FFLAGS)  briblo.f
bleblo.o: bleblo.src mg40.def cophyspar.h
	cat mg40.def bleblo.src |$(CPP) $(CPPFLAGS)  > bleblo.f
	 $(F77) $(OPTIONS) $(FFLAGS)  bleblo.f
bfrcon.o: bfrcon.src mg40.def colevel.h compi.h mgpar.h
	cat mg40.def bfrcon.src |$(CPP) $(CPPFLAGS)  > bfrcon.f
	 $(F77) $(OPTIONS) $(FFLAGS)  bfrcon.f  
bfrfix.o: bfrfix.src mg40.def
	cat mg40.def bfrfix.src |$(CPP) $(CPPFLAGS)  > bfrfix.f
	 $(F77) $(OPTIONS) $(FFLAGS)  bfrfix.f  
bfrnos.o: bfrnos.src mg40.def konsta.h
	cat mg40.def bfrnos.src |$(CPP) $(CPPFLAGS)  > bfrnos.f
	 $(F77) $(OPTIONS) $(FFLAGS)  bfrnos.f  
bfrop1.o: bfrop1.src mg40.def konsta.h
	cat mg40.def bfrop1.src |$(CPP) $(CPPFLAGS)  > bfrop1.f
	 $(F77) $(OPTIONS) $(FFLAGS)  bfrop1.f  
bfromg.o: bfromg.src mg40.def cobound.h comgrid.h mgpar.h
	cat mg40.def bfromg.src |$(CPP) $(CPPFLAGS)  > bfromg.f
	 $(F77) $(OPTIONS) $(FFLAGS)  bfromg.f 
bfrpar.o: bfrpar.src mg40.def colevel.h compi.h mgpar.h
	cat mg40.def bfrpar.src |$(CPP) $(CPPFLAGS)  > bfrpar.f
	 $(F77) $(OPTIONS) $(FFLAGS)  bfrpar.f  
bfrper.o: bfrper.src mg40.def
	cat mg40.def bfrper.src |$(CPP) $(CPPFLAGS)  > bfrper.f
	 $(F77) $(OPTIONS) $(FFLAGS)  bfrper.f  
bfrsca.o: bfrsca.src mg40.def konsta.h
	cat mg40.def bfrsca.src |$(CPP) $(CPPFLAGS)  > bfrsca.f
	 $(F77) $(OPTIONS) $(FFLAGS)  bfrsca.f 
berrhs.o: berrhs.src mg40.def
	cat mg40.def berrhs.src |$(CPP) $(CPPFLAGS)  > berrhs.f
	 $(F77) $(OPTIONS) $(FFLAGS)  berrhs.f	 	 
bl_diagnose.o: bl_diagnose.src mg40.def
	cat mg40.def bl_diagnose.src |$(CPP) $(CPPFLAGS)  > bl_diagnose.f
	 $(F77) $(OPTIONS) $(FFLAGS)  bl_diagnose.f  
blecon.o: blecon.src mg40.def colevel.h compi.h mgpar.h
	cat mg40.def blecon.src |$(CPP) $(CPPFLAGS)  > blecon.f
	 $(F77) $(OPTIONS) $(FFLAGS)  blecon.f  
blenos.o: blenos.src mg40.def konsta.h
	cat mg40.def blenos.src |$(CPP) $(CPPFLAGS)  > blenos.f
	 $(F77) $(OPTIONS) $(FFLAGS)  blenos.f  
bleop1.o: bleop1.src mg40.def
	cat mg40.def bleop1.src |$(CPP) $(CPPFLAGS)  > bleop1.f
	 $(F77) $(OPTIONS) $(FFLAGS)  bleop1.f  
bleop2.o: bleop2.src mg40.def
	cat mg40.def bleop2.src |$(CPP) $(CPPFLAGS)  > bleop2.f
	 $(F77) $(OPTIONS) $(FFLAGS)  bleop2.f  
blepar.o: blepar.src mg40.def colevel.h compi.h mgpar.h
	cat mg40.def blepar.src |$(CPP) $(CPPFLAGS)  > blepar.f
	 $(F77) $(OPTIONS) $(FFLAGS)  blepar.f  
bleper.o: bleper.src mg40.def
	cat mg40.def bleper.src |$(CPP) $(CPPFLAGS)  > bleper.f
	 $(F77) $(OPTIONS) $(FFLAGS)  bleper.f 
blesca.o: blesca.src mg40.def konsta.h
	cat mg40.def blesca.src |$(CPP) $(CPPFLAGS)  > blesca.f
	 $(F77) $(OPTIONS) $(FFLAGS)  blesca.f 		  
blesli.o: blesli.src mg40.def konsta.h
	cat mg40.def blesli.src |$(CPP) $(CPPFLAGS)  > blesli.f
	 $(F77) $(OPTIONS) $(FFLAGS)  blesli.f  
blftmg.o: blftmg.src mg40.def cobound.h comgrid.h mgpar.h
	cat mg40.def blftmg.src |$(CPP) $(CPPFLAGS)  > blftmg.f
	 $(F77) $(OPTIONS) $(FFLAGS)  blftmg.f  
bodycheck.o: bodycheck.src mg40.def cobodold.h cobody.h cogrdpro.h comgrid.h costrles.h mgpar.h
	cat mg40.def bodycheck.src |$(CPP) $(CPPFLAGS)  > bodycheck.f
	 $(F77) $(OPTIONS) $(FFLAGS)  bodycheck.f  
bodyst.o: bodyst.src mg40.def
	cat mg40.def bodyst.src |$(CPP) $(CPPFLAGS)  > bodyst.f
	 $(F77) $(OPTIONS) $(FFLAGS)  bodyst.f  
boflu.o: boflu.src mg40.def
	cat mg40.def boflu.src |$(CPP) $(CPPFLAGS)  > boflu.f
	 $(F77) $(OPTIONS) $(FFLAGS)  boflu.f  
borec.o: borec.src mg40.def
	cat mg40.def borec.src |$(CPP) $(CPPFLAGS)  > borec.f
	 $(F77) $(OPTIONS) $(FFLAGS)  borec.f  
boundmg.o: boundmg.src mg40.def
	cat mg40.def boundmg.src |$(CPP) $(CPPFLAGS)  > boundmg.f
	 $(F77) $(OPTIONS) $(FFLAGS)  boundmg.f 
boundscainter.o: boundscainter.src mg40.def
	cat mg40.def boundscainter.src |$(CPP) $(CPPFLAGS)  > boundscainter.f
	 $(F77) $(OPTIONS) $(FFLAGS)  boundscainter.f 
bparmg.o: bparmg.src mg40.def
	cat mg40.def bparmg.src |$(CPP) $(CPPFLAGS)  > bparmg.f
	 $(F77) $(OPTIONS) $(FFLAGS)  bparmg.f  
bpart.o: bpart.src
	cat mg40.def bpart.src |$(CPP) $(CPPFLAGS)  > bpart.f
	 $(F77) $(OPTIONS) $(FFLAGS)  bpart.f  
bpart_body.o: bpart_body.src
	cat mg40.def bpart_body.src |$(CPP) $(CPPFLAGS)  > bpart_body.f
	 $(F77) $(OPTIONS) $(FFLAGS)  bpart_body.f  
bphi0.o: bphi0.src mg40.def
	cat mg40.def bphi0.src |$(CPP) $(CPPFLAGS)  > bphi0.f
	 $(F77) $(OPTIONS) $(FFLAGS)  bphi0.f  
brenner.o: brenner.src 
	cat mg40.def brenner.src |$(CPP) $(CPPFLAGS)  > brenner.f
	 $(F77) $(OPTIONS) $(FFLAGS)  brenner.f  
brgtmg.o: brgtmg.src mg40.def cobound.h comgrid.h mgpar.h
	cat mg40.def brgtmg.src |$(CPP) $(CPPFLAGS)  > brgtmg.f
	 $(F77) $(OPTIONS) $(FFLAGS)  brgtmg.f  
bricon.o: bricon.src mg40.def colevel.h compi.h mgpar.h
	cat mg40.def bricon.src |$(CPP) $(CPPFLAGS)  > bricon.f
	 $(F77) $(OPTIONS) $(FFLAGS)  bricon.f  
brifix.o: brifix.src mg40.def
	cat mg40.def brifix.src |$(CPP) $(CPPFLAGS)  > brifix.f
	 $(F77) $(OPTIONS) $(FFLAGS)  brifix.f  
brinos.o: brinos.src mg40.def konsta.h
	cat mg40.def brinos.src |$(CPP) $(CPPFLAGS)  > brinos.f
	 $(F77) $(OPTIONS) $(FFLAGS)  brinos.f  
briop1.o: briop1.src mg40.def konsta.h
	cat mg40.def briop1.src |$(CPP) $(CPPFLAGS)  > briop1.f
	 $(F77) $(OPTIONS) $(FFLAGS)  briop1.f  
bripar.o: bripar.src mg40.def colevel.h compi.h mgpar.h
	cat mg40.def bripar.src |$(CPP) $(CPPFLAGS)  > bripar.f
	 $(F77) $(OPTIONS) $(FFLAGS)  bripar.f  
briper.o: briper.src mg40.def
	cat mg40.def briper.src |$(CPP) $(CPPFLAGS)  > briper.f
	 $(F77) $(OPTIONS) $(FFLAGS)  briper.f  
brisca.o: brisca.src mg40.def konsta.h
	cat mg40.def brisca.src |$(CPP) $(CPPFLAGS)  > brisca.f
	 $(F77) $(OPTIONS) $(FFLAGS)  brisca.f 			 
brisli.o: brisli.src mg40.def konsta.h
	cat mg40.def brisli.src |$(CPP) $(CPPFLAGS)  > brisli.f
	 $(F77) $(OPTIONS) $(FFLAGS)  brisli.f  
btonos.o: btonos.src mg40.def konsta.h
	cat mg40.def btonos.src |$(CPP) $(CPPFLAGS)  > btonos.f
	 $(F77) $(OPTIONS) $(FFLAGS)  btonos.f  
btoop1.o: btoop1.src mg40.def cophyspar.h konsta.h
	cat mg40.def btoop1.src |$(CPP) $(CPPFLAGS)  > btoop1.f
	 $(F77) $(OPTIONS) $(FFLAGS)  btoop1.f  
btoop2.o: btoop2.src mg40.def
	cat mg40.def btoop2.src |$(CPP) $(CPPFLAGS)  > btoop2.f
	 $(F77) $(OPTIONS) $(FFLAGS)  btoop2.f
btoop3.o: btoop3.src mg40.def
	cat mg40.def btoop3.src |$(CPP) $(CPPFLAGS)  > btoop3.f
	 $(F77) $(OPTIONS) $(FFLAGS)  btoop3.f	   
btopar.o: btopar.src mg40.def
	cat mg40.def btopar.src |$(CPP) $(CPPFLAGS)  > btopar.f
	 $(F77) $(OPTIONS) $(FFLAGS)  btopar.f  
btopmg.o: btopmg.src mg40.def cobound.h comgrid.h mgpar.h
	cat mg40.def btopmg.src |$(CPP) $(CPPFLAGS)  > btopmg.f
	 $(F77) $(OPTIONS) $(FFLAGS)  btopmg.f
btosca.o: btosca.src mg40.def konsta.h
	cat mg40.def btosca.src |$(CPP) $(CPPFLAGS)  > btosca.f
	 $(F77) $(OPTIONS) $(FFLAGS)  btosca.f  	   
btosli.o: btosli.src mg40.def konsta.h
	cat mg40.def btosli.src |$(CPP) $(CPPFLAGS)  > btosli.f
	 $(F77) $(OPTIONS) $(FFLAGS)  btosli.f  
cal_rr.o: cal_rr.src
	cat mg40.def cal_rr.src |$(CPP) $(CPPFLAGS)  > cal_rr.f
	 $(F77) $(OPTIONS) $(FFLAGS)  cal_rr.f  
cal_rrrr.o: cal_rrrr.src
	cat mg40.def cal_rrrr.src |$(CPP) $(CPPFLAGS)  > cal_rrrr.f
	 $(F77) $(OPTIONS) $(FFLAGS)  cal_rrrr.f  
calcoefd.o: calcoefd.src mg40.def
	cat mg40.def calcoefd.src |$(CPP) $(CPPFLAGS)  > calcoefd.f
	 $(F77) $(OPTIONS) $(FFLAGS)  calcoefd.f 
coefkon.o: coefkon.src mg40.def
	cat mg40.def coefkon.src |$(CPP) $(CPPFLAGS)  > coefkon.f
	 $(F77) $(OPTIONS) $(FFLAGS)  coefkon.f 
coefdiv.o: coefdiv.src mg40.def
	cat mg40.def coefdiv.src |$(CPP) $(CPPFLAGS)  > coefdiv.f
	 $(F77) $(OPTIONS) $(FFLAGS)  coefdiv.f
cotofine3d.o: cotofine3d.src mg40.def
	cat mg40.def cotofine3d.src |$(CPP) $(CPPFLAGS)  > cotofine3d.f
	 $(F77) $(OPTIONS) $(FFLAGS)  cotofine3d.f
coefscal.o: coefscal.src mg40.def
	cat mg40.def coefscal.src |$(CPP) $(CPPFLAGS)  > coefscal.f
	 $(F77) $(OPTIONS) $(FFLAGS)  coefscal.f
cop3dzero.o: cop3dzero.src mg40.def
	cat mg40.def cop3dzero.src |$(CPP) $(CPPFLAGS)  > cop3dzero.f
	 $(F77) $(OPTIONS) $(FFLAGS)  cop3dzero.f
calcoefdx.o: calcoefdx.src mg40.def
	cat mg40.def calcoefdx.src |$(CPP) $(CPPFLAGS)  > calcoefdx.f
	 $(F77) $(OPTIONS) $(FFLAGS)  calcoefdx.f  
calcoeff.o: calcoeff.src mg40.def
	cat mg40.def calcoeff.src |$(CPP) $(CPPFLAGS)  > calcoeff.f
	 $(F77) $(OPTIONS) $(FFLAGS)  calcoeff.f  
calcoeffx.o: calcoeffx.src mg40.def
	cat mg40.def calcoeffx.src |$(CPP) $(CPPFLAGS)  > calcoeffx.f
	 $(F77) $(OPTIONS) $(FFLAGS)  calcoeffx.f  
caldpdx.o: caldpdx.src mg40.def
	cat mg40.def caldpdx.src |$(CPP) $(CPPFLAGS)  > caldpdx.f
	 $(F77) $(OPTIONS) $(FFLAGS)  caldpdx.f  
calpsfak.o: calpsfak.src mg40.def
	cat mg40.def calpsfak.src |$(CPP) $(CPPFLAGS)  > calpsfak.f
	 $(F77) $(OPTIONS) $(FFLAGS)  calpsfak.f  
cap1252.o: cap1252.src mg40.def
	cat mg40.def cap1252.src |$(CPP) $(CPPFLAGS)  > cap1252.f
	 $(F77) $(OPTIONS) $(FFLAGS)  cap1252.f
cap1262.o: cap1262.src mg40.def
	cat mg40.def cap1262.src |$(CPP) $(CPPFLAGS)  > cap1262.f
	 $(F77) $(OPTIONS) $(FFLAGS)  cap1262.f
channeld.o: channeld.src mg40.def
	cat mg40.def channeld.src |$(CPP) $(CPPFLAGS)  > channeld.f
	 $(F77) $(OPTIONS) $(FFLAGS)  channeld.f  
coeffop3.o: coeffop3.src mg40.def
	cat mg40.def coeffop3.src |$(CPP) $(CPPFLAGS)  > coeffop3.f
	 $(F77) $(OPTIONS) $(FFLAGS)  coeffop3.f  
coeffop32.o: coeffop32.src mg40.def
	cat mg40.def coeffop32.src |$(CPP) $(CPPFLAGS)  > coeffop32.f
	 $(F77) $(OPTIONS) $(FFLAGS)  coeffop32.f  
compose3d.o: compose3d.src mg40.def
	cat mg40.def compose3d.src |$(CPP) $(CPPFLAGS)  > compose3d.f
	 $(F77) $(OPTIONS) $(FFLAGS)  compose3d.f  
composefield.o: composefield.src mg40.def cogrdcon.h cogrdpro.h colevel.h comgrid.h compi.h mgpar.h
	cat mg40.def composefield.src |$(CPP) $(CPPFLAGS)  > composefield.f
	 $(F77) $(OPTIONS) $(FFLAGS)  composefield.f  
composegrd.o: composegrd.src mg40.def cobody.h cobound.h cogrdcon.h cogrdpro.h colevel.h comgrid.h cophyspar.h costrles.h mgpar.h
	cat mg40.def composegrd.src |$(CPP) $(CPPFLAGS)  > composegrd.f
	 $(F77) $(OPTIONS) $(FFLAGS)  composegrd.f  
conbacmg.o: conbacmg.src mg40.def cobound.h colevel.h comgrid.h compi.h mgpar.h
	cat mg40.def conbacmg.src |$(CPP) $(CPPFLAGS)  > conbacmg.f
	 $(F77) $(OPTIONS) $(FFLAGS)  conbacmg.f  
conbapar.o: conbapar.src mg40.def colevel.h compi.h mgpar.h
	cat mg40.def conbapar.src |$(CPP) $(CPPFLAGS)  > conbapar.f
	 $(F77) $(OPTIONS) $(FFLAGS)  conbapar.f  
conbotmg.o: conbotmg.src mg40.def cobound.h colevel.h comgrid.h compi.h mgpar.h
	cat mg40.def conbotmg.src |$(CPP) $(CPPFLAGS)  > conbotmg.f
	 $(F77) $(OPTIONS) $(FFLAGS)  conbotmg.f  
confromg.o: confromg.src mg40.def cobound.h colevel.h comgrid.h compi.h mgpar.h
	cat mg40.def confromg.src |$(CPP) $(CPPFLAGS)  > confromg.f
	 $(F77) $(OPTIONS) $(FFLAGS)  confromg.f  
conlftmg.o: conlftmg.src mg40.def cobound.h colevel.h comgrid.h compi.h mgpar.h
	cat mg40.def conlftmg.src |$(CPP) $(CPPFLAGS)  > conlftmg.f
	 $(F77) $(OPTIONS) $(FFLAGS)  conlftmg.f  
connectmg.o: connectmg.src mg40.def
	cat mg40.def connectmg.src |$(CPP) $(CPPFLAGS)  > connectmg.f
	 $(F77) $(OPTIONS) $(FFLAGS)  connectmg.f  
confrpar.o: confrpar.src mg40.def colevel.h compi.h mgpar.h
	cat mg40.def confrpar.src | $(CPP) $(CPPFLAGS) > confrpar.f
	$(F77) $(OPTIONS) $(FFLAGS)  confrpar.f  
conrgtmg.o: conrgtmg.src mg40.def cobound.h colevel.h comgrid.h compi.h mgpar.h
	cat mg40.def conrgtmg.src |$(CPP) $(CPPFLAGS)  > conrgtmg.f
	 $(F77) $(OPTIONS) $(FFLAGS)  conrgtmg.f  
contopar.o: contopar.src mg40.def colevel.h compi.h mgpar.h
	cat mg40.def contopar.src |$(CPP) $(CPPFLAGS)  > contopar.f
	 $(F77) $(OPTIONS) $(FFLAGS)  contopar.f  
contopmg.o: contopmg.src mg40.def cobound.h colevel.h comgrid.h compi.h mgpar.h
	cat mg40.def contopmg.src |$(CPP) $(CPPFLAGS)  > contopmg.f
	 $(F77) $(OPTIONS) $(FFLAGS)  contopmg.f  
cop3d.o: cop3d.src mg40.def
	cat mg40.def cop3d.src |$(CPP) $(CPPFLAGS)  > cop3d.f
	 $(F77) $(OPTIONS) $(FFLAGS)  cop3d.f  
copfib.o: copfib.src mg40.def
	cat mg40.def copfib.src |$(CPP) $(CPPFLAGS)  > copfib.f
	 $(F77) $(OPTIONS) $(FFLAGS)  copfib.f  
copfic.o: copfic.src mg40.def
	cat mg40.def copfic.src |$(CPP) $(CPPFLAGS)  > copfic.f
	 $(F77) $(OPTIONS) $(FFLAGS)  copfic.f  
copyregion.o: copyregion.src mg40.def
	cat mg40.def copyregion.src |$(CPP) $(CPPFLAGS)  > copyregion.f
	 $(F77) $(OPTIONS) $(FFLAGS)  copyregion.f  
cousin.o: cousin.src mg40.def conles.h konsta.h
	cat mg40.def cousin.src |$(CPP) $(CPPFLAGS)  > cousin.f
	 $(F77) $(OPTIONS) $(FFLAGS)  cousin.f  
cubinf.o: cubinf.src mg40.def
	cat mg40.def cubinf.src |$(CPP) $(CPPFLAGS)  > cubinf.f
	 $(F77) $(OPTIONS) $(FFLAGS)  cubinf.f  
d20.o: d20.src mg40.def
	cat mg40.def d20.src |$(CPP) $(CPPFLAGS)  > d20.f
	 $(F77) $(OPTIONS) $(FFLAGS)  d20.f
deconvolv2.o: deconvolv2.src mg40.def
	cat mg40.def deconvolv2.src |$(CPP) $(CPPFLAGS)  > deconvolv2.f
	 $(F77) $(OPTIONS) $(FFLAGS)  deconvolv2.f
deconvolvsca.o: deconvolvsca.src mg40.def
	cat mg40.def deconvolvsca.src |$(CPP) $(CPPFLAGS)  > deconvolvsca.f
	 $(F77) $(OPTIONS) $(FFLAGS)  deconvolvsca.f
deibi.o: deibi.src mg40.def konsta.h
	cat mg40.def deibi.src |$(CPP) $(CPPFLAGS)  > deibi.f
	 $(F77) $(OPTIONS) $(FFLAGS)  deibi.f  
deici.o: deici.src mg40.def konsta.h
	cat mg40.def deici.src |$(CPP) $(CPPFLAGS)  > deici.f
	 $(F77) $(OPTIONS) $(FFLAGS)  deici.f  
deobi.o: deobi.src mg40.def konsta.h
	cat mg40.def deobi.src |$(CPP) $(CPPFLAGS)  > deobi.f
	 $(F77) $(OPTIONS) $(FFLAGS)  deobi.f  
deoci.o: deoci.src mg40.def konsta.h
	cat mg40.def deoci.src |$(CPP) $(CPPFLAGS)  > deoci.f
	 $(F77) $(OPTIONS) $(FFLAGS)  deoci.f  
dfdx.o: dfdx.src setsta.src mg40.def sumsta.src
	cat mg40.def dfdx.src |$(CPP) $(CPPFLAGS)  > dfdx.f
	$(F77) $(OPTIONS) $(FFLAGS)  dfdx.f
dib.o: dib.src mg40.def cogrdpro.h colevel.h compi.h mgpar.h
	cat mg40.def dib.src |$(CPP) $(CPPFLAGS)  > dib.f
	 $(F77) $(OPTIONS) $(FFLAGS)  dib.f  
dibca.o: dibca.src mg40.def cstaca.h cstadi.h cstapa.h
	cat mg40.def dibca.src |$(CPP) $(CPPFLAGS)  > dibca.f
	 $(F77) $(OPTIONS) $(FFLAGS)  dibca.f  
dibca1.o: dibca1.src mg40.def cstac1.h cstad1.h cstap1.h
	cat mg40.def dibca1.src |$(CPP) $(CPPFLAGS)  > dibca1.f
	 $(F77) $(OPTIONS) $(FFLAGS)  dibca1.f  
dibhead.o: dibhead.src mg40.def cobodold.h cogrdold.h comgrid.h costrles.h mgpar.h
	cat mg40.def dibhead.src |$(CPP) $(CPPFLAGS)  > dibhead.f
	 $(F77) $(OPTIONS) $(FFLAGS)  dibhead.f  
dic.o: dic.src mg40.def
	cat mg40.def dic.src |$(CPP) $(CPPFLAGS)  > dic.f
	 $(F77) $(OPTIONS) $(FFLAGS)  dic.f  
dichead.o: dichead.src mg40.def cobodold.h cogrdold.h comgrid.h costrles.h mgpar.h
	cat mg40.def dichead.src |$(CPP) $(CPPFLAGS)  > dichead.f
	 $(F77) $(OPTIONS) $(FFLAGS)  dichead.f  
dicrec.o: dicrec.src mg40.def corecdef.h
	cat mg40.def dicrec.src |$(CPP) $(CPPFLAGS)  > dicrec.f
	 $(F77) $(OPTIONS) $(FFLAGS)  dicrec.f  
diffcoefkobs0.o: diffcoefkobs0.src mg40.def
	cat mg40.def diffcoefkobs0.src |$(CPP) $(CPPFLAGS)  > diffcoefkobs0.f
	 $(F77) $(OPTIONS) $(FFLAGS)  diffcoefkobs0.f
digrid.o: digrid.src mg40.def clinoh.h clinou.h cobody.h cobound.h cogrdcon.h cogrdpro.h colevel.h comgrid.h compi.h cophyspar.h costrles.h mgpar.h
	cat mg40.def digrid.src |$(CPP) $(CPPFLAGS)  > digrid.f
	 $(F77) $(OPTIONS) $(FFLAGS)  digrid.f  
dimlo1.o: dimlo1.src mg40.def crefva.h cstac1.h cstad1.h cstap1.h konsta.h
	cat mg40.def dimlo1.src |$(CPP) $(CPPFLAGS)  > dimlo1.f
	 $(F77) $(OPTIONS) $(FFLAGS)  dimlo1.f  
dimlos.o: dimlos.src mg40.def crefva.h cstaca.h cstadi.h cstapa.h konsta.h
	cat mg40.def dimlos.src |$(CPP) $(CPPFLAGS)  > dimlos.f
	 $(F77) $(OPTIONS) $(FFLAGS)  dimlos.f  
dissipg.o: dissipg.src mg40.def cophyspar.h
	cat mg40.def dissipg.src |$(CPP) $(CPPFLAGS)  > dissipg.f
	 $(F77) $(OPTIONS) $(FFLAGS)  dissipg.f  
dissipga.o: dissipga.src mg40.def cophyspar.h
	cat mg40.def dissipga.src |$(CPP) $(CPPFLAGS)  > dissipga.f
	 $(F77) $(OPTIONS) $(FFLAGS)  dissipga.f
divcal.o: divcal.src mg40.def
	cat mg40.def divcal.src |$(CPP) $(CPPFLAGS)  > divcal.f
	 $(F77) $(OPTIONS) $(FFLAGS)  divcal.f  
dmixle.o: dmixle.src mg40.def conles.h konsta.h
	cat mg40.def dmixle.src |$(CPP) $(CPPFLAGS)  > dmixle.f
	 $(F77) $(OPTIONS) $(FFLAGS)  dmixle.f  
dob.o: dob.src mg40.def cogrdpro.h colevel.h compi.h konsta.h mgpar.h
	cat mg40.def dob.src |$(CPP) $(CPPFLAGS)  > dob.f
	 $(F77) $(OPTIONS) $(FFLAGS)  dob.f  
dobca.o: dobca.src mg40.def cdobca.h cstaca.h cstadi.h cstapa.h
	cat mg40.def dobca.src |$(CPP) $(CPPFLAGS)  > dobca.f
	 $(F77) $(OPTIONS) $(FFLAGS)  dobca.f  
dobca1.o: dobca1.src mg40.def cdobca.h cstac1.h cstad1.h cstap1.h
	cat mg40.def dobca1.src |$(CPP) $(CPPFLAGS)  > dobca1.f
	 $(F77) $(OPTIONS) $(FFLAGS)  dobca1.f  
dobhead.o: dobhead.src mg40.def cobody.h cobound.h cogrdcon.h cogrddef.h cogrdpro.h comgrid.h costrles.h mgpar.h
	cat mg40.def dobhead.src |$(CPP) $(CPPFLAGS)  > dobhead.f
	 $(F77) $(OPTIONS) $(FFLAGS)  dobhead.f  
dobrec.o: dobrec.src mg40.def
	cat mg40.def dobrec.src |$(CPP) $(CPPFLAGS)  > dobrec.f
	 $(F77) $(NOPTIONS) $(FFLAGS)  dobrec.f  
doc.o: doc.src mg40.def konsta.h
	cat mg40.def doc.src |$(CPP) $(CPPFLAGS)  > doc.f
	 $(F77) $(OPTIONS) $(FFLAGS)  doc.f  
dochead.o: dochead.src mg40.def cobody.h cobound.h cogrdcon.h cogrddef.h cogrdpro.h comgrid.h costrles.h mgpar.h
	cat mg40.def dochead.src |$(CPP) $(CPPFLAGS)  > dochead.f
	 $(F77) $(OPTIONS) $(FFLAGS)  dochead.f  
doenst.o: doenst.src mg40.def cdobca.h crefva.h cstaca.h cstadi.h cstapa.h konsta.h
	cat mg40.def doenst.src |$(CPP) $(CPPFLAGS)  > doenst.f
	 $(F77) $(OPTIONS) $(FFLAGS)  doenst.f  
dogrid.o: dogrid.src mg40.def clinoh.h clinou.h comgrid.h costrles.h mgpar.h
	cat mg40.def dogrid.src |$(CPP) $(CPPFLAGS)  > dogrid.f
	 $(F77) $(OPTIONS) $(FFLAGS)  dogrid.f  
dphi0.o: dphi0.src mg40.def
	cat mg40.def dphi0.src |$(CPP) $(CPPFLAGS)  > dphi0.f
	 $(F77) $(OPTIONS) $(FFLAGS)  dphi0.f  
dtmxcal.o: dtmxcal.src mg40.def
	cat mg40.def dtmxcal.src |$(CPP) $(CPPFLAGS)  > dtmxcal.f
	 $(F77) $(OPTIONS) $(FFLAGS)  dtmxcal.f  
dumpfi.o: dumpfi.src mg40.def
	cat mg40.def dumpfi.src |$(CPP) $(CPPFLAGS)  > dumpfi.f
	 $(F77) $(OPTIONS) $(FFLAGS)  dumpfi.f  
efvisc.o: efvisc.src mg40.def conles.h konsta.h
	cat mg40.def efvisc.src |$(CPP) $(CPPFLAGS)  > efvisc.f
	 $(F77) $(OPTIONS) $(FFLAGS)  efvisc.f  
enerfg.o: enerfg.src mg40.def
	cat mg40.def enerfg.src |$(CPP) $(CPPFLAGS)  > enerfg.f
	 $(F77) $(OPTIONS) $(FFLAGS)  enerfg.f  
enerfs.o: enerfs.src mg40.def conles.h konsta.h
	cat mg40.def enerfs.src |$(CPP) $(CPPFLAGS)  > enerfs.f
	 $(F77) $(OPTIONS) $(FFLAGS)  enerfs.f  
enstro.o: enstro.src mg40.def conles.h konsta.h
	cat mg40.def enstro.src |$(CPP) $(CPPFLAGS)  > enstro.f
	 $(F77) $(OPTIONS) $(FFLAGS)  enstro.f  
ernorm.o: ernorm.src mg40.def
	cat mg40.def ernorm.src |$(CPP) $(CPPFLAGS)  > ernorm.f
	 $(F77) $(OPTIONS) $(FFLAGS)  ernorm.f  
errr.o: errr.src mg40.def
	cat mg40.def errr.src |$(CPP) $(CPPFLAGS)  > errr.f
	 $(F77) $(OPTIONS) $(FFLAGS)  errr.f  
exchan.o: exchan.src mg40.def
	cat mg40.def exchan.src |$(CPP) $(CPPFLAGS)  > exchan.f
	 $(F77) $(OPTIONS) $(FFLAGS)  exchan.f  
exchansca.o: exchansca.src mg40.def
	cat mg40.def exchansca.src |$(CPP) $(CPPFLAGS)  > exchansca.f
	 $(F77) $(OPTIONS) $(FFLAGS)  exchansca.f  	 
expan1.o: expan1.src mg40.def
	cat mg40.def expan1.src |$(CPP) $(CPPFLAGS)  > expan1.f
	 $(F77) $(OPTIONS) $(FFLAGS)  expan1.f  
extendk.o: extendk.src mg40.def
	cat mg40.def extendk.src |$(CPP) $(CPPFLAGS)  > extendk.f
	$(F77) $(OPTIONS) $(FFLAGS)  extendk.f
fderfouvz.o: fderfouvz.src mg40.def
	cat mg40.def fderfouvz.src |$(CPP) $(CPPFLAGS)  > fderfouvz.f
	 $(F77) $(OPTIONS) $(FFLAGS)  fderfouvz.f 
fderfosca.o: fderfosca.src mg40.def
	cat mg40.def fderfosca.src |$(CPP) $(CPPFLAGS)  > fderfosca.f
	 $(F77) $(OPTIONS) $(FFLAGS)  fderfosca.f 
fderfouwy.o: fderfouwy.src mg40.def
	cat mg40.def fderfouwy.src |$(CPP) $(CPPFLAGS)  > fderfouwy.f
	 $(F77) $(OPTIONS) $(FFLAGS)  fderfouwy.f  
fderfouwyper.o: fderfouwyper.src mg40.def
	cat mg40.def fderfouwyper.src |$(CPP) $(CPPFLAGS)  > fderfouwyper.f
	 $(F77) $(OPTIONS) $(FFLAGS)  fderfouwyper.f  
fderfoux.o: fderfoux.src mg40.def
	cat mg40.def fderfoux.src |$(CPP) $(CPPFLAGS)  > fderfoux.f
	 $(F77) $(OPTIONS) $(FFLAGS)  fderfoux.f  
fderfouxper.o: fderfouxper.src mg40.def
	cat mg40.def fderfouxper.src |$(CPP) $(CPPFLAGS)  > fderfouxper.f
	 $(F77) $(OPTIONS) $(FFLAGS)  fderfouxper.f  
fderfovwx.o: fderfovwx.src mg40.def
	cat mg40.def fderfovwx.src |$(CPP) $(CPPFLAGS)  > fderfovwx.f
	 $(F77) $(OPTIONS) $(FFLAGS)  fderfovwx.f  
fderfovwxper.o: fderfovwxper.src mg40.def
	cat mg40.def fderfovwxper.src |$(CPP) $(CPPFLAGS)  > fderfovwxper.f
	 $(F77) $(OPTIONS) $(FFLAGS)  fderfovwxper.f  
fderfovy.o: fderfovy.src mg40.def
	cat mg40.def fderfovy.src |$(CPP) $(CPPFLAGS)  > fderfovy.f
	 $(F77) $(OPTIONS) $(FFLAGS)  fderfovy.f  
fderfovyper.o: fderfovyper.src mg40.def
	cat mg40.def fderfovyper.src |$(CPP) $(CPPFLAGS)  > fderfovyper.f
	 $(F77) $(OPTIONS) $(FFLAGS)  fderfovyper.f  
fderfowz.o: fderfowz.src mg40.def
	cat mg40.def fderfowz.src |$(CPP) $(CPPFLAGS)  > fderfowz.f
	 $(F77) $(OPTIONS) $(FFLAGS)  fderfowz.f  
fftpack.o: fftpack.src mg40.def
	cat mg40.def fftpack.src |$(CPP) $(CPPFLAGS)  > fftpack.f
	 $(F77) $(OPTIONS) $(FFLAGS)  fftpack.f  
findposition.o: findposition.src mg40.def konsta.h
	cat mg40.def findposition.src |$(CPP) $(CPPFLAGS)  > findposition.f
	 $(F77) $(OPTIONS) $(FFLAGS)  findposition.f  
filter_explizit.o: filter_explizit.src mg40.def 
	cat mg40.def filter_explizit.src |$(CPP) $(CPPFLAGS)  > filter_explizit.f
	 $(F77) $(OPTIONS) $(FFLAGS)  filter_explizit.f
filter_periodicx.o: filter_periodicx.src mg40.def 
	cat mg40.def filter_periodicx.src |$(CPP) $(CPPFLAGS)  > filter_periodicx.f
	 $(F77) $(OPTIONS) $(FFLAGS)  filter_periodicx.f
filter_periodicy.o: filter_periodicy.src mg40.def 
	cat mg40.def filter_periodicy.src |$(CPP) $(CPPFLAGS)  > filter_periodicy.f
	 $(F77) $(OPTIONS) $(FFLAGS)  filter_periodicy.f
flalfa.o: flalfa.src mg40.def
	cat mg40.def flalfa.src |$(CPP) $(CPPFLAGS)  > flalfa.f
	 $(F77) $(OPTIONS) $(FFLAGS)  flalfa.f  
flgama.o: flgama.src mg40.def
	cat mg40.def flgama.src |$(CPP) $(CPPFLAGS)  > flgama.f
	 $(F77) $(OPTIONS) $(FFLAGS)  flgama.f  
fluctuations.o: fluctuations.src mg40.def
	cat mg40.def fluctuations.src |$(CPP) $(CPPFLAGS)  > fluctuations.f
	 $(F77) $(OPTIONS) $(FFLAGS)  fluctuations.f  
fmgout.o: fmgout.src mg40.def
	cat mg40.def fmgout.src |$(CPP) $(CPPFLAGS)  > fmgout.f
	 $(F77) $(OPTIONS) $(FFLAGS)  fmgout.f  
fsensi.o: fsensi.src mg40.def
	cat mg40.def fsensi.src |$(CPP) $(CPPFLAGS)  > fsensi.f
	 $(F77) $(OPTIONS) $(FFLAGS)  fsensi.f  
fsensj.o: fsensj.src mg40.def
	cat mg40.def fsensj.src |$(CPP) $(CPPFLAGS)  > fsensj.f
	 $(F77) $(OPTIONS) $(FFLAGS)  fsensj.f  
fsensk.o: fsensk.src mg40.def
	cat mg40.def fsensk.src |$(CPP) $(CPPFLAGS)  > fsensk.f
	 $(F77) $(OPTIONS) $(FFLAGS)  fsensk.f  
funduc.o: funduc.src mg40.def conles.h conpot.h cophyspar.h mgpar.h
	cat mg40.def funduc.src |$(CPP) $(CPPFLAGS)  > funduc.f
	 $(F77) $(OPTIONS) $(FFLAGS)  funduc.f  
galirt.o: galirt.src mg40.def
	cat mg40.def galirt.src |$(CPP) $(CPPFLAGS)  > galirt.f
	 $(F77) $(OPTIONS) $(FFLAGS)  galirt.f  
getderivatives.o: getderivatives.src 
	cat mg40.def getderivatives.src |$(CPP) $(CPPFLAGS)  > getderivatives.f
	 $(F77) $(OPTIONS) $(FFLAGS)  getderivatives.f  
getvelocities.o: getvelocities.src 
	cat mg40.def getvelocities.src |$(CPP) $(CPPFLAGS)  > getvelocities.f
	 $(F77) $(OPTIONS) $(FFLAGS)  getvelocities.f  
getscalar.o: getscalar.src
	cat mg40.def getscalar.src |$(CPP) $(CPPFLAGS)  > getscalar.f
	 $(F77) $(OPTIONS) $(FFLAGS)  getscalar.f
grdbko.o: grdbko.src mg40.def
	cat mg40.def grdbko.src |$(CPP) $(CPPFLAGS)  > grdbko.f
	 $(F77) $(OPTIONS) $(FFLAGS)  grdbko.f  
grdche.o: grdche.src mg40.def konsta.h
	cat mg40.def grdche.src |$(CPP) $(CPPFLAGS)  > grdche.f
	 $(F77) $(OPTIONS) $(FFLAGS)  grdche.f  
grdchk.o: grdchk.src mg40.def cogrdcon.h cogrdpro.h comgrid.h mgpar.h
	cat mg40.def grdchk.src |$(CPP) $(CPPFLAGS)  > grdchk.f
	 $(F77) $(OPTIONS) $(FFLAGS)  grdchk.f  
grdctof.o: grdctof.src mg40.def colevel.h compi.h mgpar.h
	cat mg40.def grdctof.src |$(CPP) $(CPPFLAGS)  > grdctof.f
	 $(F77) $(OPTIONS) $(FFLAGS)  grdctof.f  
grdfit.o: grdfit.src mg40.def
	cat mg40.def grdfit.src |$(CPP) $(CPPFLAGS)  > grdfit.f
	 $(F77) $(OPTIONS) $(FFLAGS)  grdfit.f  
grdfmi.o: grdfmi.src mg40.def colevel.h compi.h mgpar.h
	cat mg40.def grdfmi.src |$(CPP) $(CPPFLAGS)  > grdfmi.f
	 $(F77) $(OPTIONS) $(FFLAGS)  grdfmi.f  
grdftoc.o: grdftoc.src mg40.def colevel.h compi.h mgpar.h
	cat mg40.def grdftoc.src |$(CPP) $(CPPFLAGS)  > grdftoc.f
	 $(F77) $(OPTIONS) $(FFLAGS)  grdftoc.f  
grdkon.o: grdkon.src mg40.def cobound.h cogrdcon.h cogrdold.h cogrdpro.h comgrid.h costrles.h mgpar.h
	cat mg40.def grdkon.src |$(CPP) $(CPPFLAGS)  > grdkon.f
	 $(F77) $(OPTIONS) $(FFLAGS)  grdkon.f  
grdkor.o: grdkor.src mg40.def konsta.h
	cat mg40.def grdkor.src |$(CPP) $(CPPFLAGS)  > grdkor.f
	 $(F77) $(OPTIONS) $(FFLAGS)  grdkor.f  
grdout.o: grdout.src mg40.def colevel.h compi.h mgpar.h
	cat mg40.def grdout.src |$(CPP) $(CPPFLAGS)  > grdout.f
	 $(F77) $(OPTIONS) $(FFLAGS)  grdout.f  
grdrec.o: grdrec.src mg40.def
	cat mg40.def grdrec.src |$(CPP) $(CPPFLAGS)  > grdrec.f
	 $(F77) $(OPTIONS) $(FFLAGS)  grdrec.f  
gridae.o: gridae.src mg40.def
	cat mg40.def gridae.src |$(CPP) $(CPPFLAGS)  > gridae.f
	 $(F77) $(OPTIONS) $(FFLAGS)  gridae.f
gsit.o: gsit.src mg40.def
	cat mg40.def gsit.src |$(CPP) $(CPPFLAGS)  > gsit.f
	 $(F77) $(OPTIONS) $(FFLAGS)  gsit.f
giteig.o: giteig.src mg40.def
	cat mg40.def giteig.src |$(CPP) $(CPPFLAGS)  > giteig.f
	 $(F77) $(OPTIONS) $(FFLAGS)  giteig.f	   
helici.o: helici.src mg40.def
	cat mg40.def helici.src |$(CPP) $(CPPFLAGS)  > helici.f
	 $(F77) $(OPTIONS) $(FFLAGS)  helici.f  
hrelom.o: hrelom.src mg40.def konsta.h
	cat mg40.def hrelom.src |$(CPP) $(CPPFLAGS)  > hrelom.f
	 $(F77) $(OPTIONS) $(FFLAGS)  hrelom.f  
htmles.o: htmles.src mg40.def
	cat mg40.def htmles.src |$(CPP) $(CPPFLAGS)  > htmles.f
	 $(F77) $(OPTIONS) $(FFLAGS)  htmles.f  
ibfield.o: ibfield.src mg40.def cobody.h cogrdpro.h conles.h cophyspar.h konsta.h mgpar.h
	cat mg40.def ibfield.src |$(CPP) $(CPPFLAGS)  > ibfield.f
	 $(F77) $(OPTIONS) $(FFLAGS)  ibfield.f  
icobody.o: icobody.src mg40.def cobound.h cogrddef.h comgrid.h mgpar.h
	cat mg40.def icobody.src |$(CPP) $(CPPFLAGS)  > icobody.f
	 $(F77) $(OPTIONS) $(FFLAGS)  icobody.f  
icobound.o: icobound.src mg40.def cobound.h cogrddef.h comgrid.h konsta.h mgpar.h
	cat mg40.def icobound.src |$(CPP) $(CPPFLAGS)  > icobound.f
	 $(F77) $(OPTIONS) $(FFLAGS)  icobound.f  
icogrdcon.o: icogrdcon.src mg40.def cogrdcon.h mgpar.h
	cat mg40.def icogrdcon.src |$(CPP) $(CPPFLAGS)  > icogrdcon.f
	 $(F77) $(OPTIONS) $(FFLAGS)  icogrdcon.f  
icogrddef.o: icogrddef.src mg40.def cobodold.h cogrdpro.h comgrid.h costrles.h mgpar.h
	cat mg40.def icogrddef.src |$(CPP) $(CPPFLAGS)  > icogrddef.f
	 $(F77) $(OPTIONS) $(FFLAGS)  icogrddef.f  
icogrdpro.o: icogrdpro.src mg40.def cobodold.h cogrdpro.h comgrid.h costrles.h mgpar.h
	cat mg40.def icogrdpro.src |$(CPP) $(CPPFLAGS)  > icogrdpro.f
	 $(F77) $(OPTIONS) $(FFLAGS)  icogrdpro.f  
icolevel.o: icolevel.src mg40.def colevel.h
	cat mg40.def icolevel.src |$(CPP) $(CPPFLAGS)  > icolevel.f
	 $(F77) $(OPTIONS) $(FFLAGS)  icolevel.f  
icomgrid.o: icomgrid.src mg40.def comgrid.h mgpar.h
	cat mg40.def icomgrid.src |$(CPP) $(CPPFLAGS)  > icomgrid.f
	 $(F77) $(OPTIONS) $(FFLAGS)  icomgrid.f  
icomgvp.o: icomgvp.src mg40.def comgrid.h comgvp.h mgpar.h
	cat mg40.def icomgvp.src |$(CPP) $(CPPFLAGS)  > icomgvp.f
	 $(F77) $(OPTIONS) $(FFLAGS)  icomgvp.f  
icophyspar.o: icophyspar.src mg40.def cophyspar.h
	cat mg40.def icophyspar.src |$(CPP) $(CPPFLAGS)  > icophyspar.f
	 $(F77) $(OPTIONS) $(FFLAGS)  icophyspar.f  
inigrid.o: inigrid.src mg40.def cobody.h cobound.h cogrdcon.h cogrdpro.h colevel.h comgrid.h compi.h cophyspar.h costrles.h mgpar.h
	cat mg40.def inigrid.src |$(CPP) $(CPPFLAGS)  > inigrid.f
	 $(F77) $(OPTIONS) $(FFLAGS)  inigrid.f  
inislice.o: inislice.src mg40.def cobody.h cogrdcon.h cogrdpro.h comgrid.h cophyspar.h mgpar.h
	cat mg40.def inislice.src |$(CPP) $(CPPFLAGS)  > inislice.f
	 $(F77) $(OPTIONS) $(FFLAGS)  inislice.f  
initsca.o: initsca.src mg40.def cobody.h cobound.h cogrdcon.h cogrdpro.h colevel.h comgrid.h compi.h cophyspar.h costrles.h mgpar.h
	cat mg40.def initsca.src |$(CPP) $(CPPFLAGS)  > initsca.f
	 $(F77) $(OPTIONS) $(FFLAGS)  initsca.f 
initorrsommer.o: initorrsommer.src mg40.def cobody.h cobound.h cogrdcon.h cogrdpro.h colevel.h comgrid.h compi.h cophyspar.h costrles.h mgpar.h
	cat mg40.def initorrsommer.src |$(CPP) $(CPPFLAGS)  > initorrsommer.f
	 $(F77) $(OPTIONS) $(FFLAGS)  initorrsommer.f 
init_part.o: init_part.src mg40.def
	cat mg40.def init_part.src |$(CPP) $(CPPFLAGS)  > init_part.f
	 $(F77) $(OPTIONS) $(FFLAGS)  init_part.f 	 	 
init_filter_periodic.o: init_filter_periodic.src mg40.def
	cat mg40.def init_filter_periodic.src |$(CPP) $(CPPFLAGS)  > init_filter_periodic.f
	 $(F77) $(OPTIONS) $(FFLAGS)  init_filter_periodic.f 	 	 
intercoef1.o: intercoef1.src mg40.def
	cat mg40.def intercoef1.src |$(CPP) $(CPPFLAGS)  > intercoef1.f
	 $(F77) $(OPTIONS) $(FFLAGS)  intercoef1.f  
intercoef2.o: intercoef2.src mg40.def
	cat mg40.def intercoef2.src |$(CPP) $(CPPFLAGS)  > intercoef2.f
	 $(F77) $(OPTIONS) $(FFLAGS)  intercoef2.f  
intercoef3.o: intercoef3.src mg40.def
	cat mg40.def intercoef3.src |$(CPP) $(CPPFLAGS)  > intercoef3.f
	 $(F77) $(OPTIONS) $(FFLAGS)  intercoef3.f  
intercoef4.o: intercoef4.src mg40.def
	cat mg40.def intercoef4.src |$(CPP) $(CPPFLAGS)  > intercoef4.f
	 $(F77) $(OPTIONS) $(FFLAGS)  intercoef4.f  
intercoef5.o: intercoef5.src mg40.def
	cat mg40.def intercoef5.src |$(CPP) $(CPPFLAGS)  > intercoef5.f
	 $(F77) $(OPTIONS) $(FFLAGS)  intercoef5.f 
intercoefkobs0.o: intercoefkobs0.src mg40.def
	cat mg40.def intercoefkobs0.src |$(CPP) $(CPPFLAGS)  > intercoefkobs0.f
	 $(F77) $(OPTIONS) $(FFLAGS)  intercoefkobs0.f 
intercoefku3.o: intercoefku3.src mg40.def
	cat mg40.def intercoefku3.src |$(CPP) $(CPPFLAGS)  > intercoefku3.f
	 $(F77) $(OPTIONS) $(FFLAGS)  intercoefku3.f  
interpolate1.o: interpolate1.src mg40.def
	cat mg40.def interpolate1.src |$(CPP) $(CPPFLAGS)  > interpolate1.f
	 $(F77) $(OPTIONS) $(FFLAGS)  interpolate1.f  
interpolate2.o: interpolate2.src mg40.def
	cat mg40.def interpolate2.src |$(CPP) $(CPPFLAGS)  > interpolate2.f
	 $(F77) $(OPTIONS) $(FFLAGS)  interpolate2.f  
interpolate3.o: interpolate3.src mg40.def
	cat mg40.def interpolate3.src |$(CPP) $(CPPFLAGS)  > interpolate3.f
	 $(F77) $(OPTIONS) $(FFLAGS)  interpolate3.f  
interpolatedi.o: interpolatedi.src mg40.def
	cat mg40.def interpolatedi.src |$(CPP) $(CPPFLAGS)  > interpolatedi.f
	 $(F77) $(OPTIONS) $(FFLAGS)  interpolatedi.f  
interpolatedj.o: interpolatedj.src mg40.def
	cat mg40.def interpolatedj.src |$(CPP) $(CPPFLAGS)  > interpolatedj.f
	 $(F77) $(OPTIONS) $(FFLAGS)  interpolatedj.f  
interpolatedk.o: interpolatedk.src mg40.def
	cat mg40.def interpolatedk.src |$(CPP) $(CPPFLAGS)  > interpolatedk.f
	 $(F77) $(OPTIONS) $(FFLAGS)  interpolatedk.f  
interpolatekj.o: interpolatekj.src mg40.def
	cat mg40.def interpolatekj.src |$(CPP) $(CPPFLAGS)  > interpolatekj.f
	 $(F77) $(OPTIONS) $(FFLAGS)  interpolatekj.f  
interpolatekk.o: interpolatekk.src mg40.def
	cat mg40.def interpolatekk.src |$(CPP) $(CPPFLAGS)  > interpolatekk.f
	 $(F77) $(OPTIONS) $(FFLAGS)  interpolatekk.f  
interpolatep.o: interpolatep.src mg40.def
	cat mg40.def interpolatep.src |$(CPP) $(CPPFLAGS)  > interpolatep.f
	$(F77) $(OPTIONS) $(FFLAGS)  interpolatep.f
interpolateuvz.o: interpolateuvz.src mg40.def
	cat mg40.def interpolateuvz.src |$(CPP) $(CPPFLAGS)  > interpolateuvz.f
	 $(F77) $(OPTIONS) $(FFLAGS)  interpolateuvz.f  
interpolateuwy.o: interpolateuwy.src mg40.def
	cat mg40.def interpolateuwy.src |$(CPP) $(CPPFLAGS)  > interpolateuwy.f
	 $(F77) $(OPTIONS) $(FFLAGS)  interpolateuwy.f
interpolatesca.o: interpolatesca.src mg40.def conles.h conpot.h cophyspar.h konsta.h cobound.h
	cat mg40.def interpolatesca.src |$(CPP) $(CPPFLAGS)  > interpolatesca.f
	 $(F77) $(OPTIONS) $(FFLAGS)  interpolatesca.f
interpolateux.o: interpolateux.src mg40.def
	cat mg40.def interpolateux.src |$(CPP) $(CPPFLAGS)  > interpolateux.f
	 $(F77) $(OPTIONS) $(FFLAGS)  interpolateux.f  
interpolateux3.o: interpolateux3.src mg40.def
	cat mg40.def interpolateux3.src |$(CPP) $(CPPFLAGS)  > interpolateux3.f
	 $(F77) $(OPTIONS) $(FFLAGS)  interpolateux3.f  
interpolatevwx.o: interpolatevwx.src mg40.def
	cat mg40.def interpolatevwx.src |$(CPP) $(CPPFLAGS)  > interpolatevwx.f
	 $(F77) $(OPTIONS) $(FFLAGS)  interpolatevwx.f  
interpolatevy.o: interpolatevy.src mg40.def
	cat mg40.def interpolatevy.src |$(CPP) $(CPPFLAGS)  > interpolatevy.f
	 $(F77) $(OPTIONS) $(FFLAGS)  interpolatevy.f  
interpolatewz.o: interpolatewz.src mg40.def
	cat mg40.def interpolatewz.src |$(CPP) $(CPPFLAGS)  > interpolatewz.f
	 $(F77) $(OPTIONS) $(FFLAGS)  interpolatewz.f  
interpolphi.o: interpolphi.src mg40.def
	cat mg40.def interpolphi.src |$(CPP) $(CPPFLAGS)  > interpolphi.f
	$(F77) $(OPTIONS) $(FFLAGS)  interpolphi.f
interuwyper.o: interuwyper.src mg40.def
	cat mg40.def interuwyper.src |$(CPP) $(CPPFLAGS)  > interuwyper.f
	 $(F77) $(OPTIONS) $(FFLAGS)  interuwyper.f  
interuxper.o: interuxper.src mg40.def
	cat mg40.def interuxper.src |$(CPP) $(CPPFLAGS)  > interuxper.f
	 $(F77) $(OPTIONS) $(FFLAGS)  interuxper.f  
intervwxper.o: intervwxper.src mg40.def
	cat mg40.def intervwxper.src |$(CPP) $(CPPFLAGS)  > intervwxper.f
	 $(F77) $(OPTIONS) $(FFLAGS)  intervwxper.f  
intervyper.o: intervyper.src mg40.def
	cat mg40.def intervyper.src |$(CPP) $(CPPFLAGS)  > intervyper.f
	 $(F77) $(OPTIONS) $(FFLAGS)  intervyper.f  
interpolate_3d_cell2.o: interpolate_3d_cell2.src
	cat mg40.def interpolate_3d_cell2.src |$(CPP) $(CPPFLAGS)  > interpolate_3d_cell2.f
	 $(F77) $(OPTIONS) $(FFLAGS)  interpolate_3d_cell2.f  
interpolate_3d_cell3.o: interpolate_3d_cell3.src
	cat mg40.def interpolate_3d_cell3.src |$(CPP) $(CPPFLAGS)  > interpolate_3d_cell3.f
	 $(F77) $(OPTIONS) $(FFLAGS)  interpolate_3d_cell3.f  
interpolate_3d_cell6.o: interpolate_3d_cell6.src
	cat mg40.def interpolate_3d_cell6.src |$(CPP) $(CPPFLAGS)  > interpolate_3d_cell6.f
	 $(F77) $(OPTIONS) $(FFLAGS)  interpolate_3d_cell6.f  
intlin.o: intlin.src mg40.def
	cat mg40.def intlin.src |$(CPP) $(CPPFLAGS)  > intlin.f
	 $(F77) $(OPTIONS) $(FFLAGS)  intlin.f  
intst.o: intst.src mg40.def
	cat mg40.def intst.src |$(CPP) $(CPPFLAGS)  > intst.f
	 $(F77) $(OPTIONS) $(FFLAGS)  intst.f  
intsta.o: intsta.src mg40.def
	cat mg40.def intsta.src |$(CPP) $(CPPFLAGS)  > intsta.f
	 $(F77) $(OPTIONS) $(FFLAGS)  intsta.f  
itinf.o: itinf.src mg40.def
	cat mg40.def itinf.src |$(CPP) $(CPPFLAGS)  > itinf.f
	 $(F77) $(OPTIONS) $(FFLAGS)  itinf.f  
itsample.o: itsample.src mg40.def cogrdcon.h cogrdpro.h colevel.h comgrid.h compi.h mgpar.h
	cat mg40.def itsample.src |$(CPP) $(CPPFLAGS)  > itsample.f
	 $(F77) $(OPTIONS) $(FFLAGS)  itsample.f  
kjidco.o: kjidco.src mg40.def
	cat mg40.def kjidco.src |$(CPP) $(CPPFLAGS)  > kjidco.f
	 $(F77) $(OPTIONS) $(FFLAGS)  kjidco.f  
kjieco.o: kjieco.src mg40.def
	cat mg40.def kjieco.src |$(CPP) $(CPPFLAGS)  > kjieco.f
	 $(F77) $(OPTIONS) $(FFLAGS)  kjieco.f  
ko2hom.o: ko2hom.src mg40.def konsta.h
	cat mg40.def ko2hom.src |$(CPP) $(CPPFLAGS)  > ko2hom.f
	 $(F77) $(OPTIONS) $(FFLAGS)  ko2hom.f  
ko2pkt.o: ko2pkt.src mg40.def konsta.h
	cat mg40.def ko2pkt.src |$(CPP) $(CPPFLAGS)  > ko2pkt.f
	 $(F77) $(OPTIONS) $(FFLAGS)  ko2pkt.f  
ko2var.o: ko2var.src mg40.def
	cat mg40.def ko2var.src |$(CPP) $(CPPFLAGS)  > ko2var.f
	 $(F77) $(OPTIONS) $(FFLAGS)  ko2var.f  
kstprt.o: kstprt.src mg40.def
	cat mg40.def kstprt.src |$(CPP) $(CPPFLAGS)  > kstprt.f
	 $(F77) $(OPTIONS) $(FFLAGS)  kstprt.f  
lesconst.o: lesconst.src mg40.def conles.h konsta.h
	cat mg40.def lesconst.src |$(CPP) $(CPPFLAGS)  > lesconst.f
	 $(F77) $(OPTIONS) $(FFLAGS)  lesconst.f  
linctl.o: linctl.src mg40.def konsta.h
	cat mg40.def linctl.src |$(CPP) $(CPPFLAGS)  > linctl.f
	 $(F77) $(OPTIONS) $(FFLAGS)  linctl.f  
lininf.o: lininf.src mg40.def cstac1.h cstad1.h cstap1.h
	cat mg40.def lininf.src |$(CPP) $(CPPFLAGS)  > lininf.f
	 $(F77) $(OPTIONS) $(FFLAGS)  lininf.f  
linout.o: linout.src mg40.def clinoh.h clinou.h
	cat mg40.def linout.src |$(CPP) $(CPPFLAGS)  > linout.f
	 $(F77) $(OPTIONS) $(FFLAGS)  linout.f  
listi6.o: listi6.src mg40.def
	cat mg40.def listi6.src |$(CPP) $(CPPFLAGS)  > listi6.f
	 $(F77) $(OPTIONS) $(FFLAGS)  listi6.f  
mamili.o: mamili.src mg40.def konsta.h
	cat mg40.def mamili.src |$(CPP) $(CPPFLAGS)  > mamili.f
	 $(F77) $(OPTIONS) $(FFLAGS)  mamili.f  
mgbasb.o: mgbasb.src mg40.def cobound.h comgrid.h mgpar.h
	cat mg40.def mgbasb.src |$(CPP) $(CPPFLAGS)  > mgbasb.f
	 $(F77) $(OPTIONS) $(FFLAGS)  mgbasb.f  
mgbasbsca.o: mgbasbsca.src mg40.def cobound.h comgrid.h mgpar.h
	cat mg40.def mgbasbsca.src |$(CPP) $(CPPFLAGS)  > mgbasbsca.f
	 $(F77) $(OPTIONS) $(FFLAGS)  mgbasbsca.f 	 
mgbftc.o: mgbftc.src mg40.def colevel.h compi.h mgpar.h
	cat mg40.def mgbftc.src |$(CPP) $(CPPFLAGS)  > mgbftc.f
	 $(F77) $(OPTIONS) $(FFLAGS)  mgbftc.f  
mgbloindset.o: mgbloindset.src mg40.def cobound.h colevel.h cophyspar.h mgpar.h
	cat mg40.def mgbloindset.src |$(CPP) $(CPPFLAGS)  > mgbloindset.f
	 $(F77) $(OPTIONS) $(FFLAGS)  mgbloindset.f  
mgboflu.o: mgboflu.src mg40.def cobound.h cogrdpro.h colevel.h compi.h cophyspar.h mgpar.h
	cat mg40.def mgboflu.src |$(CPP) $(CPPFLAGS)  > mgboflu.f
	 $(F77) $(OPTIONS) $(FFLAGS)  mgboflu.f  
mgconinf.o: mgconinf.src mg40.def cobound.h colevel.h compi.h cophyspar.h mgpar.h
	cat mg40.def mgconinf.src |$(CPP) $(CPPFLAGS)  > mgconinf.f
	 $(F77) $(OPTIONS) $(FFLAGS)  mgconinf.f  
mgcover.o: mgcover.src mg40.def konsta.h
	cat mg40.def mgcover.src |$(CPP) $(CPPFLAGS)  > mgcover.f
	 $(F77) $(OPTIONS) $(FFLAGS)  mgcover.f  
mgctof.o: mgctof.src mg40.def colevel.h compi.h mgpar.h
	cat mg40.def mgctof.src |$(CPP) $(CPPFLAGS)  > mgctof.f
	 $(F77) $(OPTIONS) $(FFLAGS)  mgctof.f  
mgdima.o: mgdima.src mg40.def comgrid.h mgpar.h
	cat mg40.def mgdima.src |$(CPP) $(CPPFLAGS)  > mgdima.f
	 $(F77) $(OPTIONS) $(FFLAGS)  mgdima.f  
mgdims.o: mgdims.src mg40.def comgrid.h mgpar.h
	cat mg40.def mgdims.src |$(CPP) $(CPPFLAGS)  > mgdims.f
	 $(F77) $(OPTIONS) $(FFLAGS)  mgdims.f  
mgdpb.o: mgdpb.src mg40.def
	cat mg40.def mgdpb.src |$(CPP) $(CPPFLAGS)  > mgdpb.f
	 $(F77) $(OPTIONS) $(FFLAGS)  mgdpb.f  
mgftoc.o: mgftoc.src mg40.def colevel.h compi.h mgpar.h
	cat mg40.def mgftoc.src |$(CPP) $(CPPFLAGS)  > mgftoc.f
	 $(F77) $(OPTIONS) $(FFLAGS)  mgftoc.f  
mggrdgen.o: mggrdgen.src mg40.def cobody.h cogrdcon.h cogrddef.h cogrdpro.h comgrid.h conles.h cophyspar.h costrles.h konsta.h mgpar.h
	cat mg40.def mggrdgen.src |$(CPP) $(CPPFLAGS)  > mggrdgen.f
	 $(F77) $(OPTIONS) $(FFLAGS)  mggrdgen.f  
mgnbrbuf.o: mgnbrbuf.src mg40.def
	cat mg40.def mgnbrbuf.src |$(CPP) $(CPPFLAGS)  > mgnbrbuf.f
	 $(F77) $(OPTIONS) $(FFLAGS)  mgnbrbuf.f  
mgnbrcheck.o: mgnbrcheck.src mg40.def cobound.h comgrid.h mgpar.h
	cat mg40.def mgnbrcheck.src |$(CPP) $(CPPFLAGS)  > mgnbrcheck.f
	 $(F77) $(OPTIONS) $(FFLAGS)  mgnbrcheck.f  
mgnbrset.o: mgnbrset.src mg40.def cobound.h colevel.h cophyspar.h mgpar.h
	cat mg40.def mgnbrset.src |$(CPP) $(CPPFLAGS)  > mgnbrset.f
	 $(F77) $(OPTIONS) $(FFLAGS)  mgnbrset.f  
mgoverlap.o: mgoverlap.src mg40.def konsta.h
	cat mg40.def mgoverlap.src |$(CPP) $(CPPFLAGS)  > mgoverlap.f
	 $(F77) $(OPTIONS) $(FFLAGS)  mgoverlap.f  
mgparcheck.o: mgparcheck.src mg40.def cobound.h cogrdcon.h comgrid.h mgpar.h
	cat mg40.def mgparcheck.src |$(CPP) $(CPPFLAGS)  > mgparcheck.f
	 $(F77) $(OPTIONS) $(FFLAGS)  mgparcheck.f  
mgparset.o: mgparset.src mg40.def cobound.h cogrdcon.h cogrdpro.h colevel.h cophyspar.h mgpar.h
	cat mg40.def mgparset.src |$(CPP) $(CPPFLAGS)  > mgparset.f
	 $(F77) $(OPTIONS) $(FFLAGS)  mgparset.f  
mgpcorr.o: mgpcorr.src mg40.def
	cat mg40.def mgpcorr.src |$(CPP) $(CPPFLAGS)  > mgpcorr.f
	 $(F77) $(OPTIONS) $(FFLAGS)  mgpcorr.f  
mgpoina.o: mgpoina.src mg40.def comgrid.h mgpar.h
	cat mg40.def mgpoina.src |$(CPP) $(CPPFLAGS)  > mgpoina.f
	 $(F77) $(OPTIONS) $(FFLAGS)  mgpoina.f  
mgpoint.o: mgpoint.src mg40.def comgrid.h mgpar.h
	cat mg40.def mgpoint.src |$(CPP) $(CPPFLAGS)  > mgpoint.f
	 $(F77) $(OPTIONS) $(FFLAGS)  mgpoint.f 
mgpoisc1.o: mgpoisc1.src mg40.def cogrdcon.h cogrdpro.h  colevel.h comgrid.h comgvp.h compi.h costrles.h mgpar.h
	cat mg40.def mgpoisc1.src |$(CPP) $(CPPFLAGS)  > mgpoisc1.f
	 $(F77) $(OPTIONS) $(FFLAGS)  mgpoisc1.f
mgpoisl1.o: mgpoisl1.src mg40.def cogrdcon.h cogrdpro.h  colevel.h comgrid.h comgvp.h compi.h costrles.h mgpar.h
	cat mg40.def mgpoisl1.src |$(CPP) $(CPPFLAGS)  > mgpoisl1.f
	 $(F77) $(OPTIONS) $(FFLAGS)  mgpoisl1.f
mgpoisit.o: mgpoisit.src mg40.def cogrdcon.h cogrdpro.h  colevel.h comgrid.h comgvp.h compi.h costrles.h mgpar.h
	cat mg40.def mgpoisit.src |$(CPP) $(CPPFLAGS)  > mgpoisit.f
	 $(F77) $(OPTIONS) $(FFLAGS)  mgpoisit.f	  
mgpsdir.o: mgpsdir.src mg40.def cogrdcon.h coksr.h colevel.h comgrid.h comgvp.h costrles.h konsta.h mgpar.h
	cat mg40.def mgpsdir.src |$(CPP) $(CPPFLAGS)  > mgpsdir.f
	 $(F77) $(OPTIONS) $(FFLAGS)  mgpsdir.f  
mgsvflu.o: mgsvflu.src mg40.def cobound.h cogrdpro.h colevel.h compi.h cophyspar.h mgpar.h
	cat mg40.def mgsvflu.src |$(CPP) $(CPPFLAGS)  > mgsvflu.f
	 $(F77) $(OPTIONS) $(FFLAGS)  mgsvflu.f  
mgtstpart.o: mgtstpart.src mg40.def copart.h 
	cat mg40.def mgtstpart.src |$(CPP) $(CPPFLAGS)  > mgtstpart.f
	 $(F77) $(NOPTIONS) $(FFLAGS)  mgtstpart.f  
mgvpc1.o: mgvpc1.src mg40.def cogrdcon.h cogrdpro.h coksr.h colevel.h comgrid.h comgvp.h compi.h costrles.h mgpar.h
	cat mg40.def mgvpc1.src |$(CPP) $(CPPFLAGS)  > mgvpc1.f
	 $(F77) $(OPTIONS) $(FFLAGS)  mgvpc1.f  
mgvpc2.o: mgvpc2.src mg40.def cogrdcon.h cogrdpro.h coksr.h colevel.h comgrid.h comgvp.h compi.h costrles.h mgpar.h
	cat mg40.def mgvpc2.src |$(CPP) $(CPPFLAGS)  > mgvpc2.f
	 $(F77) $(OPTIONS) $(FFLAGS)  mgvpc2.f  
mgvpc3.o: mgvpc3.src mg40.def cogrdcon.h cogrdpro.h coksr.h colevel.h comgrid.h comgvp.h compi.h costrles.h mgpar.h
	cat mg40.def mgvpc3.src |$(CPP) $(CPPFLAGS)  > mgvpc3.f
	 $(F77) $(OPTIONS) $(FFLAGS)  mgvpc3.f  
mgvpc4.o: mgvpc4.src mg40.def cogrdcon.h cogrdpro.h coksr.h colevel.h comgrid.h comgvp.h compi.h costrles.h mgpar.h
	cat mg40.def mgvpc4.src |$(CPP) $(CPPFLAGS)  > mgvpc4.f
	 $(F77) $(OPTIONS) $(FFLAGS)  mgvpc4.f  
mgvpc5.o: mgvpc5.src mg40.def cogrdcon.h cogrdpro.h coksr.h colevel.h comgrid.h comgvp.h compi.h costrles.h mgpar.h
	cat mg40.def mgvpc5.src |$(CPP) $(CPPFLAGS)  > mgvpc5.f
	 $(F77) $(OPTIONS) $(FFLAGS)  mgvpc5.f  
mgvpit.o: mgvpit.src mg40.def cogrdcon.h coksr.h colevel.h comgrid.h comgvp.h compi.h costrles.h mgpar.h
	cat mg40.def mgvpit.src |$(CPP) $(CPPFLAGS)  > mgvpit.f
	 $(F77) $(OPTIONS) $(FFLAGS)  mgvpit.f  
mgvpl1.o: mgvpl1.src mg40.def cogrdcon.h cogrdpro.h coksr.h colevel.h comgrid.h comgvp.h compi.h costrles.h mgpar.h
	cat mg40.def mgvpl1.src |$(CPP) $(CPPFLAGS)  > mgvpl1.f
	 $(F77) $(OPTIONS) $(FFLAGS)  mgvpl1.f  
mgvpset.o: mgvpset.src mg40.def cobody.h cobound.h cogrdcon.h cogrdpro.h colevel.h comgrid.h compi.h cophyspar.h costrles.h mgpar.h
	cat mg40.def mgvpset.src |$(CPP) $(CPPFLAGS)  > mgvpset.f
	 $(F77) $(OPTIONS) $(FFLAGS)  mgvpset.f  
mlet.o: mlet.src mg40.def clinoh.h clinou.h cobound.h cogrdcon.h cogrddef.h cogrdpro.h coksr.h colevel.h comgrid.h comgvp.h compi.h conles.h conpot.h cophyspar.h coprint.h corecdef.h costrles.h cstmg1.h cstmga.h cstsga.h konsta.h mgpar.h stat0.h stat1.h stat2.h stat3.h stat4.h stat5.h stat6.h stat7.h copart.h copart2.h
	cat mg40.def mlet.src |$(CPP) $(CPPFLAGS)  > mlet.f
	 $(F77) $(NOPTIONS) $(FFLAGS)  mlet.f  
msgvaropt.o: msgvaropt.src mg40.def cobound.h cogrdpro.h colevel.h compi.h cophyspar.h mgpar.h
	cat mg40.def msgvaropt.src |$(CPP) $(CPPFLAGS)  > msgvaropt.f
	 $(F77) $(OPTIONS) $(FFLAGS)  msgvaropt.f  
nextgrid.o: nextgrid.src mg40.def comgrid.h mgpar.h
	cat mg40.def nextgrid.src |$(CPP) $(CPPFLAGS)  > nextgrid.f
	 $(F77) $(OPTIONS) $(FFLAGS)  nextgrid.f  
omega.o: omega.src mg40.def konsta.h
	cat mg40.def omega.src |$(CPP) $(CPPFLAGS)  > omega.f
	 $(F77) $(OPTIONS) $(FFLAGS)  omega.f  
optfreq.o: optfreq.src mg40.def mgpar.h colevel.h compi.h costrles.h cobound.h
	cat mg40.def optfreq.src |$(CPP) $(CPPFLAGS)  > optfreq.f
	 $(F77) $(OPTIONS) $(FFLAGS)  optfreq.f  
output.o: output.src mg40.def
	cat mg40.def output.src |$(CPP) $(CPPFLAGS)  > output.f
	 $(F77) $(OPTIONS) $(FFLAGS)  output.f  
outres.o: outres.src mg40.def
	cat mg40.def outres.src |$(CPP) $(CPPFLAGS)  > outres.f
	 $(F77) $(OPTIONS) $(FFLAGS)  outres.f  
orrsommer.o: orrsommer.src mg40.def
	cat mg40.def orrsommer.src |$(CPP) $(CPPFLAGS)  > orrsommer.f
	 $(F77) $(OPTIONS) $(FFLAGS)  orrsommer.f  
parbacmg.o: parbacmg.src mg40.def cobound.h colevel.h comgrid.h compi.h mgpar.h
	cat mg40.def parbacmg.src |$(CPP) $(CPPFLAGS)  > parbacmg.f
	 $(F77) $(OPTIONS) $(FFLAGS)  parbacmg.f  
parbotmg.o: parbotmg.src mg40.def cobound.h colevel.h comgrid.h compi.h mgpar.h
	cat mg40.def parbotmg.src |$(CPP) $(CPPFLAGS)  > parbotmg.f
	 $(F77) $(OPTIONS) $(FFLAGS)  parbotmg.f  
parfromg.o: parfromg.src mg40.def cobound.h colevel.h comgrid.h compi.h mgpar.h
	cat mg40.def parfromg.src |$(CPP) $(CPPFLAGS)  > parfromg.f
	 $(F77) $(OPTIONS) $(FFLAGS)  parfromg.f  
parlftmg.o: parlftmg.src mg40.def cobound.h colevel.h comgrid.h compi.h mgpar.h
	cat mg40.def parlftmg.src |$(CPP) $(CPPFLAGS)  > parlftmg.f
	 $(F77) $(OPTIONS) $(FFLAGS)  parlftmg.f  
parrgtmg.o: parrgtmg.src mg40.def cobound.h colevel.h comgrid.h compi.h mgpar.h
	cat mg40.def parrgtmg.src |$(CPP) $(CPPFLAGS)  > parrgtmg.f
	 $(F77) $(OPTIONS) $(FFLAGS)  parrgtmg.f  
parset.o: parset.src mg40.def cobound.h cogrdcon.h cogrdpro.h colevel.h cophyspar.h mgpar.h
	cat mg40.def parset.src |$(CPP) $(CPPFLAGS)  > parset.f
	 $(F77) $(OPTIONS) $(FFLAGS)  parset.f  
partdiss.o: partdiss.src mg40.def
	cat mg40.def partdiss.src |$(CPP) $(CPPFLAGS)  > partdiss.f
	 $(F77) $(OPTIONS) $(FFLAGS)  partdiss.f
particle_stress.o: particle_stress.src
	cat mg40.def particle_stress.src |$(CPP) $(CPPFLAGS)  > particle_stress.f
	 $(F77) $(OPTIONS) $(FFLAGS)  particle_stress.f  
partopmg.o: partopmg.src mg40.def cobound.h colevel.h comgrid.h compi.h mgpar.h
	cat mg40.def partopmg.src |$(CPP) $(CPPFLAGS)  > partopmg.f
	 $(F77) $(OPTIONS) $(FFLAGS)  partopmg.f 
passpart.o: passpart.src mg40.def copart2.h
	cat mg40.def passpart.src |$(CPP) $(CPPFLAGS)  > passpart.f
	 $(F77) $(OPTIONS) $(FFLAGS)  passpart.f 
perturb.o: perturb.src mg40.def copart2.h
	cat mg40.def perturb.src |$(CPP) $(CPPFLAGS)  > perturb.f
	 $(F77) $(OPTIONS) $(FFLAGS)  perturb.f 
phiadd.o: phiadd.src mg40.def
	cat mg40.def phiadd.src |$(CPP) $(CPPFLAGS)  > phiadd.f
	 $(F77) $(OPTIONS) $(FFLAGS)  phiadd.f  
phididj.o: phididj.src mg40.def konsta.h
	cat mg40.def phididj.src |$(CPP) $(CPPFLAGS)  > phididj.f
	 $(F77) $(OPTIONS) $(FFLAGS)  phididj.f  
phifl2.o: phifl2.src mg40.def
	cat mg40.def phifl2.src |$(CPP) $(CPPFLAGS)  > phifl2.f
	 $(F77) $(OPTIONS) $(FFLAGS)  phifl2.f  
phifla.o: phifla.src mg40.def konsta.h
	cat mg40.def phifla.src |$(CPP) $(CPPFLAGS)  > phifla.f
	 $(F77) $(OPTIONS) $(FFLAGS)  phifla.f  
phiflu.o: phiflu.src mg40.def
	cat mg40.def phiflu.src |$(CPP) $(CPPFLAGS)  > phiflu.f
	 $(F77) $(OPTIONS) $(FFLAGS)  phiflu.f  
phimlt2.o: phimlt2.src mg40.def
	cat mg40.def phimlt2.src |$(CPP) $(CPPFLAGS) > phimlt2.f
	$(F77) $(OPTIONS) $(FFLAGS)  phimlt2.f
phimlt3.o: phimlt3.src mg40.def
	cat mg40.def phimlt3.src |$(CPP) $(CPPFLAGS) > phimlt3.f
	$(F77) $(OPTIONS) $(FFLAGS)  phimlt3.f
phimlt4.o: phimlt4.src mg40.def
	cat mg40.def phimlt4.src |$(CPP) $(CPPFLAGS) > phimlt4.f
	 $(F77) $(OPTIONS) $(FFLAGS)  phimlt4.f
phirmi.o: phirmi.src mg40.def
	cat mg40.def phirmi.src |$(CPP) $(CPPFLAGS)  > phirmi.f
	 $(F77) $(OPTIONS) $(FFLAGS)  phirmi.f  
phiske.o: phiske.src mg40.def konsta.h
	cat mg40.def phiske.src |$(CPP) $(CPPFLAGS)  > phiske.f
	 $(F77) $(OPTIONS) $(FFLAGS)  phiske.f  
plevel.o: plevel.src mg40.def
	cat mg40.def plevel.src |$(CPP) $(CPPFLAGS)  > plevel.f
	 $(F77) $(OPTIONS) $(FFLAGS)  plevel.f  
pnivea.o: pnivea.src mg40.def
	cat mg40.def pnivea.src |$(CPP) $(CPPFLAGS)  > pnivea.f
	 $(F77) $(OPTIONS) $(FFLAGS)  pnivea.f  
posrec.o: posrec.src mg40.def
	cat mg40.def posrec.src |$(CPP) $(CPPFLAGS)  > posrec.f
	 $(F77) $(OPTIONS) $(FFLAGS)  posrec.f  
pr1lin.o: pr1lin.src mg40.def
	cat mg40.def pr1lin.src |$(CPP) $(CPPFLAGS)  > pr1lin.f
	 $(F77) $(OPTIONS) $(FFLAGS)  pr1lin.f  
pr3b.o: pr3b.src mg40.def
	cat mg40.def pr3b.src |$(CPP) $(CPPFLAGS)  > pr3b.f
	 $(F77) $(OPTIONS) $(FFLAGS)  pr3b.f  
preproc_field.o: preproc_field.src mg40.def
	cat mg40.def preproc_field.src |$(CPP) $(CPPFLAGS)  > preproc_field.f
	 $(F77) $(OPTIONS) $(FFLAGS)  preproc_field.f  
preproc_field_scalar.o: preproc_field_scalar.src mg40.def
	cat mg40.def preproc_field_scalar.src |$(CPP) $(CPPFLAGS)  > preproc_field_scalar.f
	 $(F77) $(OPTIONS) $(FFLAGS)  preproc_field_scalar.f  
preproc_field_wcomp.o: preproc_field_wcomp.src mg40.def
	cat mg40.def preproc_field_wcomp.src |$(CPP) $(CPPFLAGS)  > preproc_field_wcomp.f
	 $(F77) $(OPTIONS) $(FFLAGS)  preproc_field_wcomp.f  
preproc_pressure.o: preproc_pressure.src mg40.def
	cat mg40.def preproc_pressure.src |$(CPP) $(CPPFLAGS)  > preproc_pressure.f
	 $(F77) $(OPTIONS) $(FFLAGS)  preproc_pressure.f  
pr3v.o: pr3v.src mg40.def
	cat mg40.def pr3v.src |$(CPP) $(CPPFLAGS)  > pr3v.f
	 $(F77) $(OPTIONS) $(FFLAGS)  pr3v.f  
printe.o: printe.src mg40.def
	cat mg40.def printe.src |$(CPP) $(CPPFLAGS)  > printe.f
	 $(F77) $(OPTIONS) $(FFLAGS)  printe.f  
printvel.o: printvel.src mg40.def
	cat mg40.def printvel.src |$(CPP) $(CPPFLAGS)  > printvel.f
	 $(F77) $(OPTIONS) $(FFLAGS)  printvel.f  
prle3e.o: prle3e.src mg40.def crefvc.h crefvd.h cstaca.h cstadi.h cstapa.h
	cat mg40.def prle3e.src |$(CPP) $(CPPFLAGS)  > prle3e.f
	 $(F77) $(OPTIONS) $(FFLAGS)  prle3e.f  
prle3l.o: prle3l.src mg40.def crefvc.h crefvd.h cstac1.h cstad1.h cstap1.h
	cat mg40.def prle3l.src |$(CPP) $(CPPFLAGS)  > prle3l.f
	 $(F77) $(OPTIONS) $(FFLAGS)  prle3l.f  
prle3m.o: prle3m.src mg40.def coprint.h
	cat mg40.def prle3m.src |$(CPP) $(CPPFLAGS)  > prle3m.f
	 $(F77) $(OPTIONS) $(FFLAGS)  prle3m.f  
prolong1.o: prolong1.src mg40.def
	cat mg40.def prolong1.src |$(CPP) $(CPPFLAGS)  > prolong1.f
	 $(F77) $(OPTIONS) $(FFLAGS)  prolong1.f  
prolong2.o: prolong2.src mg40.def
	cat mg40.def prolong2.src |$(CPP) $(CPPFLAGS)  > prolong2.f
	 $(F77) $(OPTIONS) $(FFLAGS)  prolong2.f  
random.o: random.src mg40.def
	cat mg40.def random.src |$(CPP) $(CPPFLAGS)  > random.f
	 $(F77) $(NOPTIONS) $(FFLAGS)  random.f  
random_functions.o: random_functions.src mg40.def
	cat mg40.def random_functions.src |$(CPP) $(CPPFLAGS)  > random_functions.f
	$(F77) $(OPTIONS) $(FFLAGS)  random_functions.f
random_numbers.o: random_numbers.src mg40.def
	cat mg40.def random_numbers.src |$(CPP) $(CPPFLAGS)  > random_numbers.f
	$(F77) $(OPTIONS) $(FFLAGS)  random_numbers.f
randwert.o: randwert.src mg40.def
	cat mg40.def randwert.src |$(CPP) $(CPPFLAGS)  > randwert.f
	 $(F77) $(OPTIONS) $(FFLAGS)  randwert.f  
randwertan.o: randwertan.src mg40.def
	cat mg40.def randwertan.src |$(CPP) $(CPPFLAGS)  > randwertan.f
	 $(F77) $(OPTIONS) $(FFLAGS)  randwertan.f  
randwertand.o: randwertand.src mg40.def
	cat mg40.def randwertand.src |$(CPP) $(CPPFLAGS)  > randwertand.f
	 $(F77) $(OPTIONS) $(FFLAGS)  randwertand.f  
randwertd.o: randwertd.src mg40.def
	cat mg40.def randwertd.src |$(CPP) $(CPPFLAGS)  > randwertd.f
	 $(F77) $(OPTIONS) $(FFLAGS)  randwertd.f  
randwertdw.o: randwertdw.src mg40.def
	cat mg40.def randwertdw.src |$(CPP) $(CPPFLAGS)  > randwertdw.f
	 $(F77) $(OPTIONS) $(FFLAGS)  randwertdw.f  
randwertij.o: randwertij.src mg40.def
	cat mg40.def randwertij.src |$(CPP) $(CPPFLAGS)  > randwertij.f
	 $(F77) $(OPTIONS) $(FFLAGS)  randwertij.f  
randwertijan.o: randwertijan.src mg40.def
	cat mg40.def randwertijan.src |$(CPP) $(CPPFLAGS)  > randwertijan.f
	 $(F77) $(OPTIONS) $(FFLAGS)  randwertijan.f  
randwertijand.o: randwertijand.src mg40.def
	cat mg40.def randwertijand.src |$(CPP) $(CPPFLAGS)  > randwertijand.f
	 $(F77) $(OPTIONS) $(FFLAGS)  randwertijand.f  
randwertijandw.o: randwertijandw.src mg40.def
	cat mg40.def randwertijandw.src |$(CPP) $(CPPFLAGS)  > randwertijandw.f
	 $(F77) $(OPTIONS) $(FFLAGS)  randwertijandw.f  
randwertijd.o: randwertijd.src mg40.def
	cat mg40.def randwertijd.src |$(CPP) $(CPPFLAGS)  > randwertijd.f
	 $(F77) $(OPTIONS) $(FFLAGS)  randwertijd.f  
randwertijdw.o: randwertijdw.src mg40.def
	cat mg40.def randwertijdw.src |$(CPP) $(CPPFLAGS)  > randwertijdw.f
	 $(F77) $(OPTIONS) $(FFLAGS)  randwertijdw.f  
randwertijw.o: randwertijw.src mg40.def
	cat mg40.def randwertijw.src |$(CPP) $(CPPFLAGS)  > randwertijw.f
	 $(F77) $(OPTIONS) $(FFLAGS)  randwertijw.f  
randwertw.o: randwertw.src mg40.def
	cat mg40.def randwertw.src |$(CPP) $(CPPFLAGS)  > randwertw.f
	 $(F77) $(OPTIONS) $(FFLAGS)  randwertw.f  
read_particles.o: read_particles.src 
	cat mg40.def read_particles.src |$(CPP) $(CPPFLAGS)  > read_particles.f
	 $(F77) $(OPTIONS) $(FFLAGS)  read_particles.f 
read_particles_ini.o: read_particles_ini.src
	cat mg40.def read_particles_ini.src |$(CPP) $(CPPFLAGS)  > read_particles_ini.f
	 $(F77) $(OPTIONS) $(FFLAGS)  read_particles_ini.f 
readsliced.o: readsliced.src mg40.def cogrdcon.h cogrdpro.h colevel.h compi.h mgpar.h
	cat mg40.def readsliced.src |$(CPP) $(CPPFLAGS)  > readsliced.f
	 $(F77) $(OPTIONS) $(FFLAGS)  readsliced.f  
recout.o: recout.src mg40.def
	cat mg40.def recout.src |$(CPP) $(CPPFLAGS)  > recout.f
	 $(F77) $(OPTIONS) $(FFLAGS)  recout.f  
rezipd.o: rezipd.src mg40.def
	cat mg40.def rezipd.src |$(CPP) $(CPPFLAGS)  > rezipd.f
	 $(F77) $(OPTIONS) $(FFLAGS)  rezipd.f  
rkjiae.o: rkjiae.src mg40.def
	cat mg40.def rkjiae.src |$(CPP) $(CPPFLAGS)  > rkjiae.f
	 $(F77) $(OPTIONS) $(FFLAGS)  rkjiae.f  
rmcomment.o: rmcomment.src mg40.def
	cat mg40.def rmcomment.src |$(CPP) $(CPPFLAGS)  > rmcomment.f
	 $(F77) $(OPTIONS) $(FFLAGS)  rmcomment.f  
sbbc12.o: sbbc12.src mg40.def colevel.h compi.h conles.h konsta.h mgpar.h
	cat mg40.def sbbc12.src |$(CPP) $(CPPFLAGS)  > sbbc12.f
	 $(F77) $(OPTIONS) $(FFLAGS)  sbbc12.f  
sbcoun.o: sbcoun.src mg40.def conles.h konsta.h
	cat mg40.def sbcoun.src |$(CPP) $(CPPFLAGS)  > sbcoun.f
	 $(F77) $(OPTIONS) $(FFLAGS)  sbcoun.f  
sbcub.o: sbcub.src mg40.def conles.h konsta.h
	cat mg40.def sbcub.src |$(CPP) $(CPPFLAGS)  > sbcub.f
	 $(F77) $(OPTIONS) $(FFLAGS)  sbcub.f  
sbhemi.o: sbhemi.src mg40.def conles.h konsta.h
	cat mg40.def sbhemi.src |$(CPP) $(CPPFLAGS)  > sbhemi.f
	 $(F77) $(OPTIONS) $(FFLAGS)  sbhemi.f  
sbzyl.o: sbzyl.src mg40.def conles.h konsta.h
	cat mg40.def sbzyl.src |$(CPP) $(CPPFLAGS)  > sbzyl.f
	 $(F77) $(OPTIONS) $(FFLAGS)  sbzyl.f  
searchindex.o: searchindex.src mg40.def
	cat mg40.def searchindex.src |$(CPP) $(CPPFLAGS)  > searchindex.f
	 $(F77) $(OPTIONS) $(FFLAGS)  searchindex.f  
sel0.o: sel0.src mg40.def
	cat mg40.def sel0.src |$(CPP) $(CPPFLAGS)  > sel0.f
	 $(F77) $(OPTIONS) $(FFLAGS)  sel0.f  
sel1.o: sel1.src mg40.def
	cat mg40.def sel1.src |$(CPP) $(CPPFLAGS)  > sel1.f
	 $(F77) $(OPTIONS) $(FFLAGS)  sel1.f  
sel2.o: sel2.src mg40.def
	cat mg40.def sel2.src |$(CPP) $(CPPFLAGS)  > sel2.f
	 $(F77) $(OPTIONS) $(FFLAGS)  sel2.f  
sel3.o: sel3.src mg40.def
	cat mg40.def sel3.src |$(CPP) $(CPPFLAGS)  > sel3.f
	 $(F77) $(OPTIONS) $(FFLAGS)  sel3.f  
sel4.o: sel4.src mg40.def
	cat mg40.def sel4.src |$(CPP) $(CPPFLAGS)  > sel4.f
	 $(F77) $(OPTIONS) $(FFLAGS)  sel4.f  
sel5.o: sel5.src mg40.def
	cat mg40.def sel5.src |$(CPP) $(CPPFLAGS)  > sel5.f
	 $(F77) $(OPTIONS) $(FFLAGS)  sel5.f  
sel6.o: sel6.src mg40.def
	cat mg40.def sel6.src |$(CPP) $(CPPFLAGS)  > sel6.f
	 $(F77) $(OPTIONS) $(FFLAGS)  sel6.f  
sel7.o: sel7.src mg40.def
	cat mg40.def sel7.src |$(CPP) $(CPPFLAGS)  > sel7.f
	 $(F77) $(OPTIONS) $(FFLAGS)  sel7.f  
sel8.o: sel8.src mg40.def
	cat mg40.def sel8.src |$(CPP) $(CPPFLAGS)  > sel8.f
	 $(F77) $(OPTIONS) $(FFLAGS)  sel8.f  
sel11.o: sel11.src mg40.def
	cat mg40.def sel11.src |$(CPP) $(CPPFLAGS)  > sel11.f
	 $(F77) $(OPTIONS) $(FFLAGS)  sel11.f  
sel12.o: sel12.src mg40.def
	cat mg40.def sel12.src |$(CPP) $(CPPFLAGS)  > sel12.f
	 $(F77) $(OPTIONS) $(FFLAGS)  sel12.f  
sel13.o: sel13.src mg40.def
	cat mg40.def sel13.src |$(CPP) $(CPPFLAGS)  > sel13.f
	 $(F77) $(OPTIONS) $(FFLAGS)  sel13.f  
sel14.o: sel14.src mg40.def
	cat mg40.def sel14.src |$(CPP) $(CPPFLAGS)  > sel14.f
	 $(F77) $(OPTIONS) $(FFLAGS)  sel14.f  
sel15.o: sel15.src mg40.def
	cat mg40.def sel15.src |$(CPP) $(CPPFLAGS)  > sel15.f
	 $(F77) $(OPTIONS) $(FFLAGS)  sel15.f
sel17.o: sel17.src mg40.def
	cat mg40.def sel17.src |$(CPP) $(CPPFLAGS)  > sel17.f
	 $(F77) $(OPTIONS) $(FFLAGS)  sel17.f
sel20.o: sel20.src mg40.def
	cat mg40.def sel20.src |$(CPP) $(CPPFLAGS)  > sel20.f
	 $(F77) $(OPTIONS) $(FFLAGS)  sel20.f
selaco.o: selaco.src mg40.def colevel.h compi.h mgpar.h
	cat mg40.def selaco.src |$(CPP) $(CPPFLAGS)  > selaco.f
	 $(F77) $(OPTIONS) $(FFLAGS)  selaco.f  
selauf.o: selauf.src mg40.def
	cat mg40.def selauf.src |$(CPP) $(CPPFLAGS)  > selauf.f
	 $(F77) $(OPTIONS) $(FFLAGS)  selauf.f  
seleca.o: seleca.src mg40.def
	cat mg40.def seleca.src |$(CPP) $(CPPFLAGS)  > seleca.f
	 $(F77) $(OPTIONS) $(FFLAGS)  seleca.f  
seleci.o: seleci.src mg40.def
	cat mg40.def seleci.src |$(CPP) $(CPPFLAGS)  > seleci.f
	 $(F77) $(OPTIONS) $(FFLAGS)  seleci.f  
setbackold.o: setbackold.src mg40.def conles.h cophyspar.h konsta.h mgpar.h
	cat mg40.def setbackold.src |$(CPP) $(CPPFLAGS)  > setbackold.f
	 $(F77) $(OPTIONS) $(FFLAGS)  setbackold.f  
setbloindex.o: setbloindex.src mg40.def konsta.h
	cat mg40.def setbloindex.src |$(CPP) $(CPPFLAGS)  > setbloindex.f
	 $(F77) $(OPTIONS) $(FFLAGS)  setbloindex.f  
setcobone.o: setcobone.src mg40.def cobound.h cogrdcon.h comgrid.h mgpar.h
	cat mg40.def setcobone.src |$(CPP) $(CPPFLAGS)  > setcobone.f
	 $(F77) $(OPTIONS) $(FFLAGS)  setcobone.f  
setcobound.o: setcobound.src mg40.def cobound.h cogrdcon.h comgrid.h mgpar.h
	cat mg40.def setcobound.src |$(CPP) $(CPPFLAGS)  > setcobound.f
	 $(F77) $(OPTIONS) $(FFLAGS)  setcobound.f  
setcolevel.o: setcolevel.src mg40.def colevel.h
	cat mg40.def setcolevel.src |$(CPP) $(CPPFLAGS)  > setcolevel.f
	 $(F77) $(OPTIONS) $(FFLAGS)  setcolevel.f  
setcomgrid.o: setcomgrid.src mg40.def cogrdpro.h colevel.h comgrid.h compi.h mgpar.h
	cat mg40.def setcomgrid.src |$(CPP) $(CPPFLAGS)  > setcomgrid.f
	 $(F77) $(OPTIONS) $(FFLAGS)  setcomgrid.f  
setgeovp.o: setgeovp.src mg40.def
	cat mg40.def setgeovp.src |$(CPP) $(CPPFLAGS)  > setgeovp.f
	 $(F77) $(OPTIONS) $(FFLAGS)  setgeovp.f  
setgrd.o: setgrd.src mg40.def
	cat mg40.def setgrd.src |$(CPP) $(CPPFLAGS)  > setgrd.f
	 $(F77) $(OPTIONS) $(FFLAGS)  setgrd.f  
setgrdpro.o: setgrdpro.src mg40.def cogrddef.h cogrdpro.h mgpar.h
	cat mg40.def setgrdpro.src |$(CPP) $(CPPFLAGS)  > setgrdpro.f
	 $(F77) $(OPTIONS) $(FFLAGS)  setgrdpro.f  
setid8.o: setid8.src mg40.def
	cat mg40.def setid8.src |$(CPP) $(CPPFLAGS)  > setid8.f
	 $(F77) $(OPTIONS) $(FFLAGS)  setid8.f  
setidx.o: setidx.src mg40.def
	cat mg40.def setidx.src |$(CPP) $(CPPFLAGS)  > setidx.f
	 $(F77) $(OPTIONS) $(FFLAGS)  setidx.f  
setkon.o: setkon.src mg40.def colevel.h compi.h konsta.h mgpar.h
	cat mg40.def setkon.src |$(CPP) $(CPPFLAGS)  > setkon.f
	 $(F77) $(OPTIONS) $(FFLAGS)  setkon.f  
setmpi.o: setmpi.src mg40.def cogrdpro.h colevel.h comgrid.h compi.h mgpar.h
	cat mg40.def setmpi.src |$(CPP) $(CPPFLAGS)  > setmpi.f
	 $(F77) $(OPTIONS) $(FFLAGS)  setmpi.f  
setrandom.o: setrandom.src mg40.def
	cat mg40.def setrandom.src |$(CPP) $(CPPFLAGS)  > setrandom.f
	 $(F77) $(OPTIONS) $(FFLAGS)  setrandom.f  
setref.o: setref.src mg40.def
	cat mg40.def setref.src |$(CPP) $(CPPFLAGS)  > setref.f
	 $(F77) $(OPTIONS) $(FFLAGS)  setref.f  
sets.o: sets.src mg40.def konsta.h
	cat mg40.def sets.src |$(CPP) $(CPPFLAGS)  > sets.f
	 $(F77) $(OPTIONS) $(FFLAGS)  sets.f  
setslice.o: setslice.src mg40.def cobound.h cogrdcon.h cogrdpro.h comgrid.h cophyspar.h mgpar.h
	cat mg40.def setslice.src |$(CPP) $(CPPFLAGS)  > setslice.f
	 $(F77) $(OPTIONS) $(FFLAGS)  setslice.f  
setst1.o: setst1.src mg40.def cstac1.h cstad1.h cstap1.h
	cat mg40.def setst1.src |$(CPP) $(CPPFLAGS)  > setst1.f
	 $(F77) $(OPTIONS) $(FFLAGS)  setst1.f  
setsta.o: setsta.src mg40.def cstaca.h cstadi.h cstapa.h
	cat mg40.def setsta.src |$(CPP) $(CPPFLAGS)  > setsta.f
	 $(F77) $(OPTIONS) $(FFLAGS)  setsta.f 	  
setthreads.o: setthreads.src mg40.def cogrdpro.h coksr.h colevel.h mgpar.h
	cat mg40.def setthreads.src |$(CPP) $(CPPFLAGS)  > setthreads.f
	 $(F77) $(OPTIONS) $(FFLAGS)  setthreads.f 
sipiter.o: sipiter.src mg40.def
	cat mg40.def sipiter.src | $(CPP) $(CPPFLAGS)  > sipiter.f
	 $(F77) $(OPTIONS) $(FFLAGS)  sipiter.f
siplu.o: siplu.src mg40.def
	cat mg40.def siplu.src |$(CPP) $(CPPFLAGS)  > siplu.f
	 $(F77) $(OPTIONS) $(FFLAGS)  siplu.f	  
skonle.o: skonle.src mg40.def conles.h conpot.h
	cat mg40.def skonle.src |$(CPP) $(CPPFLAGS)  > skonle.f
	 $(F77) $(OPTIONS) $(FFLAGS)  skonle.f  
slice1d.o: slice1d.src mg40.def
	cat mg40.def slice1d.src |$(CPP) $(CPPFLAGS)  > slice1d.f
	 $(F77) $(OPTIONS) $(FFLAGS)  slice1d.f  
slice3d.o: slice3d.src mg40.def
	cat mg40.def slice3d.src |$(CPP) $(CPPFLAGS)  > slice3d.f
	 $(F77) $(OPTIONS) $(FFLAGS)  slice3d.f  
slicefield.o: slicefield.src mg40.def cogrdcon.h cogrdpro.h colevel.h comgrid.h compi.h mgpar.h
	cat mg40.def slicefield.src |$(CPP) $(CPPFLAGS)  > slicefield.f
	 $(F77) $(OPTIONS) $(FFLAGS)  slicefield.f  
slicegeo.o: slicegeo.src mg40.def cobody.h cogrdcon.h cogrdpro.h colevel.h comgrid.h compi.h cophyspar.h mgpar.h
	cat mg40.def slicegeo.src |$(CPP) $(CPPFLAGS)  > slicegeo.f
	 $(F77) $(OPTIONS) $(FFLAGS)  slicegeo.f  
slicegrd.o: slicegrd.src mg40.def cobody.h cobound.h cogrdcon.h cogrdpro.h colevel.h comgrid.h cophyspar.h costrles.h mgpar.h
	cat mg40.def slicegrd.src |$(CPP) $(CPPFLAGS)  > slicegrd.f
	 $(F77) $(OPTIONS) $(FFLAGS)  slicegrd.f  
smooth_stress.o: smooth_stress.src
	cat mg40.def smooth_stress.src |$(CPP) $(CPPFLAGS)  > smooth_stress.f
	 $(F77) $(OPTIONS) $(FFLAGS)  smooth_stress.f  
srmsgs.o: srmsgs.src mg40.def
	cat mg40.def srmsgs.src |$(CPP) $(CPPFLAGS)  > srmsgs.f
	 $(F77) $(OPTIONS) $(FFLAGS)  srmsgs.f  
stabed.o: stabed.src mg40.def
	cat mg40.def stabed.src |$(CPP) $(CPPFLAGS)  > stabed.f
	 $(F77) $(OPTIONS) $(FFLAGS)  stabed.f  
startcontrol.o: startcontrol.src mg40.def
	cat mg40.def startcontrol.src |$(CPP) $(CPPFLAGS)  > startcontrol.f
	 $(F77) $(OPTIONS) $(FFLAGS)  startcontrol.f  
stmim1.o: stmim1.src mg40.def cstac1.h cstad1.h cstap1.h
	cat mg40.def stmim1.src |$(CPP) $(CPPFLAGS)  > stmim1.f
	 $(F77) $(OPTIONS) $(FFLAGS)  stmim1.f  
stmimp.o: stmimp.src mg40.def cstaca.h cstadi.h cstapa.h
	cat mg40.def stmimp.src |$(CPP) $(CPPFLAGS)  > stmimp.f
	 $(F77) $(OPTIONS) $(FFLAGS)  stmimp.f  
stmnew.o: stmnew.src mg40.def
	cat mg40.def stmnew.src |$(CPP) $(CPPFLAGS)  > stmnew.f
	 $(F77) $(OPTIONS) $(FFLAGS)  stmnew.f  
stmnli.o: stmnli.src mg40.def
	cat mg40.def stmnli.src |$(CPP) $(CPPFLAGS)  > stmnli.f
	 $(F77) $(OPTIONS) $(FFLAGS)  stmnli.f  
stmnrm.o: stmnrm.src mg40.def
	cat mg40.def stmnrm.src |$(CPP) $(CPPFLAGS)  > stmnrm.f
	 $(F77) $(OPTIONS) $(FFLAGS)  stmnrm.f  
strcheck.o: strcheck.src mg40.def cobodold.h cogrddef.h cogrdold.h cogrdpro.h comgrid.h costrles.h mgpar.h
	cat mg40.def strcheck.src |$(CPP) $(CPPFLAGS)  > strcheck.f
	 $(F77) $(OPTIONS) $(FFLAGS)  strcheck.f  
stress.o: stress.src
	cat mg40.def stress.src |$(CPP) $(CPPFLAGS)  > stress.f
	 $(F77) $(OPTIONS) $(FFLAGS)  stress.f  
strles.o: strles.src mg40.def clinoh.h clinou.h cobody.h cobound.h cogrdcon.h cogrddef.h cogrdpro.h colevel.h comgrid.h compi.h cophyspar.h coprint.h corecdef.h costrles.h konsta.h mgpar.h copart.h copart2.h
	cat mg40.def strles.src |$(CPP) $(CPPFLAGS)  > strles.f
	 $(F77) $(OPTIONS) $(FFLAGS)  strles.f  
strset.o: strset.src mg40.def cobodold.h cogrddef.h cogrdold.h cogrdpro.h comgrid.h costrles.h mgpar.h
	cat mg40.def strset.src |$(CPP) $(CPPFLAGS)  > strset.f
	 $(F77) $(OPTIONS) $(FFLAGS)  strset.f  
sum2ar.o: sum2ar.src mg40.def
	cat mg40.def sum2ar.src |$(CPP) $(CPPFLAGS)  > sum2ar.f
	 $(F77) $(OPTIONS) $(FFLAGS)  sum2ar.f  
sumaan.o: sumaan.src mg40.def cstaca.h cstadi.h cstapa.h
	cat mg40.def sumaan.src |$(CPP) $(CPPFLAGS)  > sumaan.f
	 $(F77) $(OPTIONS) $(FFLAGS)  sumaan.f  
sumst1.o: sumst1.src mg40.def cstac1.h cstad1.h cstap1.h
	cat mg40.def sumst1.src |$(CPP) $(CPPFLAGS)  > sumst1.f
	 $(F77) $(OPTIONS) $(FFLAGS)  sumst1.f  
sumsta.o: sumsta.src mg40.def cstaca.h cstadi.h cstapa.h konsta.h
	cat mg40.def sumsta.src |$(CPP) $(CPPFLAGS)  > sumsta.f
	 $(F77) $(OPTIONS) $(FFLAGS)  sumsta.f  
sv3d.o: sv3d.src mg40.def conles.h cophyspar.h konsta.h mgpar.h
	cat mg40.def sv3d.src |$(CPP) $(CPPFLAGS)  > sv3d.f
	 $(F77) $(OPTIONS) $(FFLAGS)  sv3d.f  
sveipr.o: sveipr.src mg40.def conles.h cophyspar.h konsta.h mgpar.h
	cat mg40.def sveipr.src |$(CPP) $(CPPFLAGS)  > sveipr.f
	 $(F77) $(OPTIONS) $(FFLAGS)  sveipr.f  
sveiprsca.o: sveiprsca.src mg40.def conles.h cophyspar.h konsta.h mgpar.h
	cat mg40.def sveiprsca.src |$(CPP) $(CPPFLAGS)  > sveiprsca.f
	 $(F77) $(OPTIONS) $(FFLAGS)  sveiprsca.f  	 
svflu.o: svflu.src mg40.def
	cat mg40.def svflu.src |$(CPP) $(CPPFLAGS)  > svflu.f
	 $(F77) $(OPTIONS) $(FFLAGS)  svflu.f  
svle1.o: svle1.src mg40.def conles.h cophyspar.h konsta.h mgpar.h
	cat mg40.def svle1.src |$(CPP) $(CPPFLAGS)  > svle1.f
	 $(F77) $(NOPTIONS) $(FFLAGS)  svle1.f 
svle1sca.o: svle1sca.src mg40.def conles.h cophyspar.h konsta.h mgpar.h
	cat mg40.def svle1sca.src |$(CPP) $(CPPFLAGS)  > svle1sca.f
	 $(F77) $(OPTIONS) $(FFLAGS)  svle1sca.f  	  
svrec.o: svrec.src mg40.def
	cat mg40.def svrec.src |$(CPP) $(CPPFLAGS)  > svrec.f
	 $(F77) $(OPTIONS) $(FFLAGS)  svrec.f  
swcle3d.o: swcle3d.src mg40.def conles.h conpot.h cophyspar.h konsta.h mgpar.h
	cat mg40.def swcle3d.src |$(CPP) $(CPPFLAGS)  > swcle3d.f
	 $(F77) $(OPTIONS) $(FFLAGS)  swcle3d.f  
swclesca.o: swclesca.src mg40.def conles.h conpot.h cophyspar.h konsta.h mgpar.h
	cat mg40.def swclesca.src |$(CPP) $(CPPFLAGS)  > swclesca.f
	 $(F77) $(OPTIONS) $(FFLAGS)  swclesca.f 	 
taufg.o: taufg.src mg40.def
	cat mg40.def taufg.src |$(CPP) $(CPPFLAGS)  > taufg.f
	 $(F77) $(OPTIONS) $(FFLAGS)  taufg.f  
taufm.o: taufm.src mg40.def konsta.h
	cat mg40.def taufm.src |$(CPP) $(CPPFLAGS)  > taufm.f
	 $(F77) $(OPTIONS) $(FFLAGS)  taufm.f  
taufs.o: taufs.src mg40.def konsta.h
	cat mg40.def taufs.src |$(CPP) $(CPPFLAGS)  > taufs.f
	$(F77) $(OPTIONS) $(FFLAGS)  taufs.f
taunn.o: taunn.src mg40.def 
	cat mg40.def taunn.src |$(CPP) $(CPPFLAGS)  > taunn.f
	 $(F77) $(OPTIONS) $(FFLAGS)  taunn.f  
thomasi.o: thomasi.src mg40.def
	cat mg40.def thomasi.src |$(CPP) $(CPPFLAGS)  > thomasi.f
	 $(F77) $(OPTIONS) $(FFLAGS)  thomasi.f  
thomasj.o: thomasj.src mg40.def
	cat mg40.def thomasj.src |$(CPP) $(CPPFLAGS)  > thomasj.f
	 $(F77) $(OPTIONS) $(FFLAGS)  thomasj.f  
thomask.o: thomask.src mg40.def
	cat mg40.def thomask.src |$(CPP) $(CPPFLAGS)  > thomask.f
	 $(F77) $(OPTIONS) $(FFLAGS)  thomask.f  
trizyk.o: trizyk.src mg40.def
	cat mg40.def trizyk.src |$(CPP) $(CPPFLAGS)  > trizyk.f
	 $(F77) $(OPTIONS) $(FFLAGS)  trizyk.f  
trizyk3d.o: trizyk3d.src mg40.def
	cat mg40.def trizyk3d.src |$(CPP) $(CPPFLAGS)  > trizyk3d.f
	 $(F77) $(OPTIONS) $(FFLAGS)  trizyk3d.f
trizyk3dj.o: trizyk3dj.src mg40.def
	cat mg40.def trizyk3dj.src |$(CPP) $(CPPFLAGS)  > trizyk3dj.f
	 $(F77) $(OPTIONS) $(FFLAGS)  trizyk3dj.f  
tst1g.o: tst1g.src mg40.def clinoh.h clinou.h cobound.h cogrdcon.h cogrdpro.h coksr.h colevel.h comgrid.h comgvp.h compi.h cophyspar.h coprint.h costrles.h mgpar.h copart.h
	cat mg40.def tst1g.src |$(CPP) $(CPPFLAGS)  > tst1g.f
	 $(F77) $(OPTIONS) $(FFLAGS)  tst1g.f  
tst3rk.o: tst3rk.src mg40.def clinou.h cobound.h cogrdcon.h cogrdpro.h colevel.h comgrid.h comgvp.h compi.h cophyspar.h coprint.h costrles.h mgpar.h copart.h copart2.h
	cat mg40.def tst3rk.src |$(CPP) $(CPPFLAGS)  > tst3rk.f
	 $(F77) $(OPTIONS) $(FFLAGS)  tst3rk.f  	 
tstle2.o: tstle2.src mg40.def cophyspar.h konsta.h
	cat mg40.def tstle2.src |$(CPP) $(CPPFLAGS)  > tstle2.f
	 $(F77) $(OPTIONS) $(FFLAGS)  tstle2.f  
tstle4.o: tstle4.src mg40.def cophyspar.h konsta.h
	cat mg40.def tstle4.src |$(CPP) $(CPPFLAGS)  > tstle4.f
	 $(F77) $(OPTIONS) $(FFLAGS)  tstle4.f 	 
tstorient.o: tstorient.src 
	cat mg40.def tstorient.src |$(CPP) $(CPPFLAGS)  > tstorient.f
	 $(F77) $(NOPTIONS) $(FFLAGS)  tstorient.f  
tstpos.o: tstpos.src 
	cat mg40.def tstpos.src |$(CPP) $(CPPFLAGS)  > tstpos.f
	 $(F77) $(OPTIONS) $(FFLAGS)  tstpos.f  
tstsca.o: tstsca.src mg40.def cophyspar.h konsta.h mgpar.h cobound.h
	cat mg40.def tstsca.src |$(CPP) $(CPPFLAGS)  > tstsca.f
	 $(F77) $(OPTIONS) $(FFLAGS)  tstsca.f  	
tstsca4.o: tstsca4.src mg40.def cophyspar.h konsta.h mgpar.h cobound.h
	cat mg40.def tstsca4.src |$(CPP) $(CPPFLAGS)  > tstsca4.f
	 $(F77) $(OPTIONS) $(FFLAGS)  tstsca4.f 	  
ubulk.o: ubulk.src mg40.def
	cat mg40.def ubulk.src |$(CPP) $(CPPFLAGS)  > ubulk.f
	 $(F77) $(OPTIONS) $(FFLAGS)  ubulk.f  
vplesc.o: vplesc.src mg40.def
	cat mg40.def vplesc.src |$(CPP) $(CPPFLAGS)  > vplesc.f
	 $(F77) $(OPTIONS) $(FFLAGS)  vplesc.f  
vprbc.o: vprbc.src mg40.def
	cat mg40.def vprbc.src |$(CPP) $(CPPFLAGS)  > vprbc.f
	 $(F77) $(OPTIONS) $(FFLAGS)  vprbc.f  
vshift.o: vshift.src mg40.def konsta.h
	cat mg40.def vshift.src |$(CPP) $(CPPFLAGS)  > vshift.f
	 $(F77) $(OPTIONS) $(FFLAGS)  vshift.f 
wallkobdiffs0a.o: wallkobdiffs0a.src mg40.def
	cat mg40.def wallkobdiffs0a.src |$(CPP) $(CPPFLAGS)  > wallkobdiffs0a.f
	 $(F77) $(OPTIONS) $(FFLAGS)  wallkobdiffs0a.f
wallkobdiffs0e.o: wallkobdiffs0e.src mg40.def
	cat mg40.def wallkobdiffs0e.src |$(CPP) $(CPPFLAGS)  > wallkobdiffs0e.f
	 $(F77) $(OPTIONS) $(FFLAGS)  wallkobdiffs0e.f 
wifak.o: wifak.src mg40.def
	cat mg40.def wifak.src |$(CPP) $(CPPFLAGS)  > wifak.f
	 $(F77) $(OPTIONS) $(FFLAGS)  wifak.f 
wifaksca.o: wifaksca.src mg40.def
	cat mg40.def wifaksca.src |$(CPP) $(CPPFLAGS)  > wifaksca.f
	 $(F77) $(OPTIONS) $(FFLAGS)  wifaksca.f	  
wnsint.o: wnsint.src mg40.def conles.h konsta.h
	cat mg40.def wnsint.src |$(CPP) $(CPPFLAGS)  > wnsint.f
	 $(F77) $(OPTIONS) $(FFLAGS)  wnsint.f  
wrigeo.o: wrigeo.src mg40.def
	cat mg40.def wrigeo.src |$(CPP) $(CPPFLAGS)  > wrigeo.f
	 $(F77) $(OPTIONS) $(FFLAGS)  wrigeo.f  
wrirec.o: wrirec.src mg40.def
	cat mg40.def wrirec.src |$(CPP) $(CPPFLAGS)  > wrirec.f
	 $(F77) $(OPTIONS) $(FFLAGS)  wrirec.f  
write_dupart.o: write_dupart.src 
	cat mg40.def write_dupart.src |$(CPP) $(CPPFLAGS)  > write_dupart.f
	 $(F77) $(OPTIONS) $(FFLAGS)  write_dupart.f 
write3d.o: write3d.src
	cat mg40.def write3d.src |$(CPP) $(CPPFLAGS)  > write3d.f
	 $(F77) $(OPTIONS) $(FFLAGS)  write3d.f 
write_particles.o: write_particles.src 
	cat mg40.def write_particles.src |$(CPP) $(CPPFLAGS)  > write_particles.f
	 $(F77) $(OPTIONS) $(FFLAGS)  write_particles.f  
write_particles_fin.o: write_particles_fin.src
	cat mg40.def write_particles_fin.src |$(CPP) $(CPPFLAGS)  > write_particles_fin.f
	 $(F77) $(OPTIONS) $(FFLAGS)  write_particles_fin.f
write_passpart.o: write_passpart.src
	cat mg40.def write_passpart.src |$(CPP) $(CPPFLAGS)  > write_passpart.f
	 $(F77) $(OPTIONS) $(FFLAGS)  write_passpart.f
writesliced.o: writesliced.src mg40.def cogrdcon.h cogrdpro.h colevel.h compi.h mgpar.h
	cat mg40.def writesliced.src |$(CPP) $(CPPFLAGS)  > writesliced.f
	 $(F77) $(OPTIONS) $(FFLAGS)  writesliced.f  
wrivaropt.o: wrivaropt.src mg40.def cobound.h mgpar.h
	cat mg40.def wrivaropt.src |$(CPP) $(CPPFLAGS)  > wrivaropt.f
	 $(F77) $(OPTIONS) $(FFLAGS)  wrivaropt.f  
wssint.o: wssint.src mg40.def conles.h konsta.h
	cat mg40.def wssint.src |$(CPP) $(CPPFLAGS)  > wssint.f
	 $(F77) $(OPTIONS) $(FFLAGS)  wssint.f  
xtractxy.o: xtractxy.src mg40.def
	cat mg40.def xtractxy.src |$(CPP) $(CPPFLAGS)  > xtractxy.f
	 $(F77) $(OPTIONS) $(FFLAGS)  xtractxy.f  
xtractyz.o: xtractyz.src mg40.def
	cat mg40.def xtractyz.src |$(CPP) $(CPPFLAGS)  > xtractyz.f
	 $(F77) $(OPTIONS) $(FFLAGS)  xtractyz.f  
xtractyzf.o: xtractyzf.src mg40.def
	cat mg40.def xtractyzf.src |$(CPP) $(CPPFLAGS)  > xtractyzf.f
	 $(F77) $(OPTIONS) $(FFLAGS)  xtractyzf.f  
zbrent.o: zbrent.src mg40.def
	cat mg40.def zbrent.src |$(CPP) $(CPPFLAGS)  > zbrent.f
	 $(F77) $(OPTIONS) $(FFLAGS)  zbrent.f  
dyn.o: dyn.src mg40.def
	cat mg40.def dyn.src |$(CPP) $(CPPFLAGS)  > dyn.f
	$(F77) $(OPTIONS) $(FFLAGS)  dyn.f 
inc.o: inc.src mg40.def
	cat mg40.def inc.src |$(CPP) $(CPPFLAGS)  > inc.f
	$(F77) $(OPTIONS) $(FFLAGS)  inc.f
inch.o: inch.src mg40.def
	cat mg40.def inch.src |$(CPP) $(CPPFLAGS)  > inch.f
	$(F77) $(OPTIONS) $(FFLAGS)  inch.f
leon.o: leon.src mg40.def
	cat mg40.def leon.src |$(CPP) $(CPPFLAGS)  > leon.f
	$(F77) $(OPTIONS) $(FFLAGS)  leon.f
model.o: model.src mg40.def
	cat mg40.def model.src |$(CPP) $(CPPFLAGS)  > model.f
	$(F77) $(OPTIONS) $(FFLAGS)  model.f
tauinc.o: tauinc.src mg40.def
	cat mg40.def tauinc.src |$(CPP) $(CPPFLAGS)  > tauinc.f
	$(F77) $(OPTIONS) $(FFLAGS)  tauinc.f
conrgtmgbp.o: conrgtmgbp.src mg40.def
	cat mg40.def conrgtmgbp.src |$(CPP) $(CPPFLAGS)  > conrgtmgbp.f
	$(F77) $(OPTIONS) $(FFLAGS)  conrgtmgbp.f
confromgbp.o: confromgbp.src mg40.def
	cat mg40.def confromgbp.src |$(CPP) $(CPPFLAGS)  > confromgbp.f
	$(F77) $(OPTIONS) $(FFLAGS)  confromgbp.f
conbacmgbp.o: conbacmgbp.src mg40.def
	cat mg40.def conbacmgbp.src |$(CPP) $(CPPFLAGS)  > conbacmgbp.f
	$(F77) $(OPTIONS) $(FFLAGS)  conbacmgbp.f
conlftmgbp.o: conlftmgbp.src mg40.def
	cat mg40.def conlftmgbp.src |$(CPP) $(CPPFLAGS)  > conlftmgbp.f
	$(F77) $(OPTIONS) $(FFLAGS)  conlftmgbp.f
briconbp.o: briconbp.src mg40.def
	cat mg40.def briconbp.src |$(CPP) $(CPPFLAGS)  > briconbp.f
	$(F77) $(OPTIONS) $(FFLAGS)  briconbp.f
bleconbp.o: bleconbp.src mg40.def
	cat mg40.def bleconbp.src |$(CPP) $(CPPFLAGS)  > bleconbp.f
	$(F77) $(OPTIONS) $(FFLAGS)  bleconbp.f
bfrconbp.o: bfrconbp.src mg40.def
	cat mg40.def bfrconbp.src |$(CPP) $(CPPFLAGS)  > bfrconbp.f
	$(F77) $(OPTIONS) $(FFLAGS)  bfrconbp.f
bbaconbp.o: bbaconbp.src mg40.def
	cat mg40.def bbaconbp.src |$(CPP) $(CPPFLAGS)  > bbaconbp.f
	$(F77) $(OPTIONS) $(FFLAGS)  bbaconbp.f
inigridgeo.o: inigridgeo.src mg40.def cobody.h cobound.h cogrdcon.h cogrdpro.h colevel.h comgrid.h compi.h cophyspar.h costrles.h mgpar.h
	cat mg40.def inigridgeo.src |$(CPP) $(CPPFLAGS)  > inigridgeo.f
	$(F77) $(OPTIONS) $(FFLAGS)  inigridgeo.f
inigridvar.o: inigridvar.src mg40.def cobody.h cobound.h cogrdcon.h cogrdpro.h colevel.h comgrid.h compi.h cophyspar.h costrles.h mgpar.h
	cat mg40.def inigridvar.src |$(CPP) $(CPPFLAGS)  > inigridvar.f
	$(F77) $(OPTIONS) $(FFLAGS)  inigridvar.f
blockbp.o: blockbp.src mg40.def mgpar.h comgrid.h colevel.h cogrdcon.h compi.h cogrdpro.h cobody.h
	cat mg40.def blockbp.src |$(CPP) $(CPPFLAGS)  > blockbp.f
	$(F77) $(OPTIONS) $(FFLAGS)  blockbp.f
printdata.o: printdata.src mg40.def
	cat mg40.def printdata.src |$(CPP) $(CPPFLAGS) > printdata.f
	$(F77) $(OPTIONS) $(FFLAGS) printdata.f
#radb2.o:
#	$(F77) $(OPTIONS) $(FFLAGS) radb2.f
#radb3.o:
#	$(F77) $(OPTIONS) $(FFLAGS) radb3.f
#radb4.o:
#	$(F77) $(OPTIONS) $(FFLAGS) radb4.f
#radb5.o:
#	$(F77) $(OPTIONS) $(FFLAGS) radb5.f
#radbg.o:
#	$(F77) $(OPTIONS) $(FFLAGS) radbg.f
#radf2.o:
#	$(F77) $(OPTIONS) $(FFLAGS) radf2.f
#radf3.o:
#	$(F77) $(OPTIONS) $(FFLAGS) radf3.f
#radf4.o:
#	$(F77) $(OPTIONS) $(FFLAGS) radf4.f
#radf5.o:
#	$(F77) $(OPTIONS) $(FFLAGS) radf5.f
#radfg.o:
#	$(F77) $(OPTIONS) $(FFLAGS) radfg.f
#rfftb.o:
#	$(F77) $(OPTIONS) $(FFLAGS) rfftb.f
#rfftb1.o:
#	$(F77) $(OPTIONS) $(FFLAGS) rfftb1.f
#rfftf.o:
#	$(F77) $(OPTIONS) $(FFLAGS) rfftf.f
#rfftf1.o:
#	$(F77) $(OPTIONS) $(FFLAGS) rfftf1.f
#rffti.o:
#	$(F77) $(OPTIONS) $(FFLAGS) rffti.f
#rffti1.o:
#	$(F77) $(OPTIONS) $(FFLAGS) rffti1.f


clean : 
	rm -f *.o *.f

clobber : clean
	rm -f mg40.exe *~

