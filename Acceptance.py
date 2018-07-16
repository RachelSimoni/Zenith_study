#!/usr/bin/python

import math
import ROOT

rootSafe = []

def DrawResiduals(th_1,th_2,label):
    c4 = ROOT.TCanvas(str('Residuals_') + label,label,1000,500)
    rootSafe.append(c4)
    
    residuals12,projection12 = Residuals(th_1,th_2,True)
    rootSafe.append(residuals12)
    rootSafe.append(projection12)
    
    c4.Divide(2,1)
    c4.cd(1)
    residuals12.SetLineColor(8)
    residuals12.SetMarkerStyle(20)
    residuals12.SetMarkerColor(8)
    residuals12.SetStats(False)
    residuals12.Draw('p')
    residuals12.GetYaxis().SetRangeUser(-10.2,10.2)
    ROOT.gPad.SetTickx()
    ROOT.gPad.SetTicky()    
    ROOT.gPad.SetGridx()
    ROOT.gPad.SetGridy()    
    ROOT.gPad.Update()
    
    c4.cd(2)
    projection12.SetLineColor(8)
    projection12.SetLineWidth(3)
    projection12.SetFillColorAlpha(8,0.2)
    projection12.SetStats(False)
    projection12.Draw()
    ROOT.gPad.SetTickx()
    ROOT.gPad.SetTicky()    
    ROOT.gPad.SetGridx()
    ROOT.gPad.SetGridy()    
    ROOT.gPad.Update()
    c4.cd(1)
    lat = ROOT.TLatex()
    rootSafe.append(lat)
    lat.SetTextSize(0.08)
    lat.DrawLatex(0.2,4.,label)
    ROOT.gPad.Update()

def Residuals(hist1,hist2,percent=False):
    clone = hist2.Clone()
    clone.Reset()
    newName = str(hist1.GetName())
    newName += str(hist2.GetName())
    projection = ROOT.TH1F(newName+'_proj',';#sigma;Entries',10,-5,5)
    if percent:
        projection.Delete()
        projection = ROOT.TH1F(newName+'_proj',';Percent;Entries',20,-20,20)
    if not percent:
        clone.GetYaxis().SetTitle('#sigma')
    else:
        clone.GetYaxis().SetTitle('Percent')
    for bin in xrange(1,hist1.GetNbinsX()+1):
        cont1 = hist1.GetBinContent(bin)
        err1 = hist1.GetBinError(bin)
        cont2 = hist2.GetBinContent(bin)
        err2 = hist2.GetBinError(bin)
        diff = cont2-cont1
        err = math.sqrt(err1*err1 + err2*err2)
        if not percent:
            clone.SetBinContent(bin,diff/err)
            projection.Fill(diff/err)
        else:
            average = (cont1+cont2)/2.
            clone.SetBinContent(bin,100.*diff/average)
            projection.Fill(100.*diff/average)
    return clone,projection

def DrawLine(angle):
    xend = 5. * math.cos(angle*math.pi/180.)
    yend = 5. * math.sin(angle*math.pi/180.)
    line = ROOT.TLine(0.,0.,xend,yend)
    line.SetLineWidth(3)
    line.SetLineStyle(10)
    line.Draw()
    rootSafe.append(line)

h1 = ROOT.TH1F('hRA',';#Delta RA ( deg );Events',100,-5.,5.)
h2 = ROOT.TH1F('hDec',';#Delta Decl ( deg );Events',100,-5.,5.)
hAcc = ROOT.TH2F('acc',';#Delta RA ( deg );#Delta Decl ( deg )',500,-5.,5.,500,-5.,5.)

hZen = ROOT.TProfile2D('avgZen',';#Delta RA ( deg );#Delta Decl ( deg );#Delta Zenith ( deg )',100,-5.,5.,100,-5.,5.)
hZen.SetStats(False)
hExclusion = ROOT.TH2F('ex',';#Delta RA ( deg );#Delta Decl ( deg )',100,-5.,5.,100,-5.,5.)

th_1 = ROOT.TH1F('th_1',';r^{2} ( deg^{2} );Events',36,0.,9.)
th_2 = ROOT.TH1F('th_2',';r^{2} ( deg^{2} );Events',36,0.,9.)
th_3 = ROOT.TH1F('th_3',';r^{2} ( deg^{2} );Events',36,0.,9.)
th_4 = ROOT.TH1F('th_4',';r^{2} ( deg^{2} );Events',36,0.,9.)
th_5 = ROOT.TH1F('th_5',';r^{2} ( deg^{2} );Events',36,0.,9.)
th_6 = ROOT.TH1F('th_6',';r^{2} ( deg^{2} );Events',36,0.,9.)

energies = ROOT.TH1F('energies','Energies;E ( TeV );Events',1000,0.,100.)
zenith_angles = ROOT.TH1F('zenith','Zenith angles;Angle ( deg );Events',900,0.,90.)

area_out = ROOT.TH1F('area_out','Out area;r ( deg );Area',10,0.,5.)
area_in = ROOT.TH1F('area_in','In area;r ( deg );Area',10,0.,5.)

for line in open('InputFiles/EventDump_all_85_4ter_W4.txt','r'):
#for line in open('InputFiles/EventDump_SafeOnly_85_4ter_W4.txt','r'):
    info = line.split()
    if '#' in str(info[0]): continue
    obsZen = float(info[6])
    obsRA = float(info[9])
    obsDecl = float(info[10])
    cosDecl = math.cos(obsDecl * math.pi / 180.)
    ra = float(info[0])
    dec = float(info[1])
    isFB = int(info[2])
    isPKS = int(info[3])
    energy = float(info[4])
    zenith = float(info[5])
    dRA = (ra-obsRA) * cosDecl
    dDec = dec-obsDecl
    dist = math.sqrt(dRA*dRA + dDec*dDec)
    h1.Fill(dRA)
    h2.Fill(dDec)
    angle = 180./math.pi * math.atan2(dDec,dRA) # -20 -80 -140. 160. 100. 40. -20.
    if angle < -20. and angle > -80.:
        th_1.Fill(dist*dist)
    elif angle < -80. and angle > -140.:
        th_2.Fill(dist*dist)
    elif angle < -140. or angle > 160.:
        th_3.Fill(dist*dist)
    elif angle < 160. and angle > 100.:
        th_4.Fill(dist*dist)
    elif angle < 100. and angle > 40.:
        th_5.Fill(dist*dist)
    elif angle < 40. and angle > -20.:
        th_6.Fill(dist*dist)
        
    energies.Fill(energy)
    zenith_angles.Fill(zenith)
    # try quick cut on energy and see if gradient goes away
    # if energy < 0.7:
    #     continue
    if zenith > 30.:
        continue
#    if not isFB:
    hAcc.Fill(dRA,dDec)
    if not isFB:
        hZen.Fill(dRA,dDec,obsZen-zenith)
    if isPKS:
        hExclusion.Fill(dRA,dDec)

# for xbin in xrange(1,hExclusion.GetNbinsX()+1):
#     ra = hExclusion.GetXaxis().GetBinCenter(xbin)
#     for ybin in xrange(1,hExclusion.GetNbinsY()+1):
#         content = hExclusion.GetBinContent(xbin,ybin)
#         dec = hExclusion.GetYaxis().GetBinCenter(ybin)
#         angle = 180./math.pi * math.atan2(dec,ra)
#         dist = math.sqrt(ra*ra + dec*dec)
#         if angle > -68. and angle < 150.:
#             if not (content > 0.):
#                 area_out.Fill(dist)
#         else:
#             if not (content > 0.):
#                 area_in.Fill(dist)

c2 = ROOT.TCanvas('c2','c2',700,700)

c3_1 = ROOT.TCanvas('c3_1','c3',700,700)
c3_2 = ROOT.TCanvas('c3_2','c3',700,700)

c1 = ROOT.TCanvas('c1','Acceptance 1',900,600)
c1.Divide(3,2)
c1.cd(1)
h2.SetLineColor(2)
h2.Draw()
h1.Draw('same')
ROOT.gPad.Update()

c2.cd()
hAcc.SetMaximum(40.)
hAcc.SetStats(False)
hAcc.Draw('colz')

lat = ROOT.TLatex()
lat.SetTextSize(0.08)
lat.DrawLatex(3.,-3,'#color[2]{1}')
lat.DrawLatex(-1.5,-4,'#color[2]{2}')
lat.DrawLatex(-4.,-0.5,'#color[2]{3}')
lat.DrawLatex(-3.,3.,'#color[2]{4}')
lat.DrawLatex(1.,4.,'#color[2]{5}')
lat.DrawLatex(3.,0.5,'#color[2]{6}')


circ1 = ROOT.TEllipse(0.,0.,1.0)
circ1.SetFillStyle(0)
circ1.Draw()
circ2 = ROOT.TEllipse(0.,0.,1.5)
circ2.SetFillStyle(0)
circ2.Draw()
circ3 = ROOT.TEllipse(0.,0.,2.0)
circ3.SetFillStyle(0)
circ3.Draw()
circ4 = ROOT.TEllipse(0.,0.,2.5)
circ4.SetFillStyle(0)
circ4.Draw()

# DrawLine(80.)
# DrawLine(150.)

DrawLine(40.)
DrawLine(100.)
DrawLine(160.) 
DrawLine(-20.)
DrawLine(-80.)
DrawLine(-140.)
ROOT.gPad.Update()


# residuals12,projection12 = Residuals(th_1,th_2,False)
# residuals42,projection42 = Residuals(th_4,th_2,False)
# residuals52,projection52 = Residuals(th_5,th_2,False)
# residuals62,projection62 = Residuals(th_6,th_2,False)

# residuals13,projection13 = Residuals(th_1,th_3,False)
# residuals43,projection43 = Residuals(th_4,th_3,False)
# residuals53,projection53 = Residuals(th_5,th_3,False)
# residuals63,projection63 = Residuals(th_6,th_3,False)

DrawResiduals(th_1,th_2,'2 wrt. 1')
DrawResiduals(th_4,th_2,'2 wrt. 4')
DrawResiduals(th_5,th_2,'2 wrt. 5')
DrawResiduals(th_6,th_2,'2 wrt. 6')

DrawResiduals(th_1,th_3,'3 wrt. 1')
DrawResiduals(th_4,th_3,'3 wrt. 4')
DrawResiduals(th_5,th_3,'3 wrt. 5')
DrawResiduals(th_6,th_3,'3 wrt. 6')

th_6.SetLineWidth(2)
th_5.SetLineWidth(2)
th_4.SetLineWidth(2)
th_3.SetLineWidth(2)
th_2.SetLineWidth(2)
th_1.SetLineWidth(2)
th_5.SetLineColor(8)
th_4.SetLineColor(2)
th_3.SetLineColor(9)
th_2.SetLineColor(1)
th_1.SetLineColor(ROOT.kOrange)

c3_1.cd()
leg1 = ROOT.TLegend(0.5,0.4,0.85,0.85)
leg1.SetFillStyle(1001)
leg1.AddEntry(th_3,'Bubble, 3','l')
leg1.AddEntry(th_1,'1','l')
leg1.AddEntry(th_4,'4','l')
leg1.AddEntry(th_5,'5','l')
leg1.AddEntry(th_6,'6','l')

th_2.SetStats(False)
th_3.SetStats(False)

th_3.Draw('e hist')
th_6.Draw('same e hist')
th_5.Draw('same e hist')
th_4.Draw('same e hist')
th_1.Draw('same e hist')
leg1.Draw('same')

ROOT.gPad.SetTickx()
ROOT.gPad.SetTicky()

ROOT.gPad.Update()

c3_2.cd()
th_2.Draw('e hist')
# th_6.Draw('same e hist')
th_5.Draw('same e hist')
# th_4.Draw('same e hist')
# th_1.Draw('same e hist')

leg2 = ROOT.TLegend(0.5,0.4,0.85,0.85)
leg2.SetFillStyle(1001)
leg2.AddEntry(th_2,'Bubble, 2','l')
leg2.AddEntry(th_1,'1','l')
leg2.AddEntry(th_4,'4','l')
leg2.AddEntry(th_5,'5','l')
leg2.AddEntry(th_6,'6','l')
leg2.Draw('same')
ROOT.gPad.SetTickx()
ROOT.gPad.SetTicky()

ROOT.gPad.Update()

c1.cd(5)
hExclusion.Draw('colz')

c1.cd(3)
zenith_angles.Draw()
ROOT.gPad.SetLogy()
c1.cd(4)
energies.Draw()
ROOT.gPad.SetLogy()
c1.cd(6)
ROOT.gPad.SetRightMargin(0.2)
ROOT.gPad.SetTopMargin(0.2)
ROOT.gPad.SetBottomMargin(0.2)
ROOT.gPad.SetLeftMargin(0.2)

hZen.Draw('colz')
hZen.GetZaxis().SetRangeUser(-3.9,3.9)
ROOT.gPad.Update()
