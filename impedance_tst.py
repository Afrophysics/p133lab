import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as opt
import scipy.stats as stat

#raw data of impedance vs frequency and phase
fils = ['boxf.txt', 'boxf_1.txt', 'boxc.txt','boxa.txt']
#fils = ['boxc.txt']
#setup parameters in a list to run a loop of tests

#the theoretical models you want to test
fntns = ['p[0] + ((-1)**(0.5))*(p[1]*x - (1/(x*p[2])))', 
         '1/((1/p[0]) - ((-1)**(0.5))*(1/(x*p[1] - (1/(x*p[2])))))',
         '(p[0]**2 + (x*p[1])**2)/(p[0] + ((-1)**(0.5))*(x*(p[2]*(p[0]**2 + (x*p[1])**2)) - p[1]))',
         '(p[0]**2 + (1/(x*p[2]))**2)/(p[0] + ((-1)**(0.5))*((1/(x*p[2])) - ((p[0]**2 + (1/(x*p[2]))**2)/(x*p[1]))))',
        '((1/p[0]) - ((-1)**(0.5))*(x*p[2] - (1/(x*p[1]))))/((1/p[0])**2 + (p[2]*x - (1/(x*p[1])))**2)',
         '+ 1/((-1)**(0.5) *x*p[2])',
         '+ ((-1)**(0.5) *x*p[1])',
         '+ p[0]',
         '+ p[0]+ ((-1)**(0.5) *x*p[1])+ 1/((-1)**(0.5) *x*p[2])',
         '1/(+ 1/(p[0])+ 1/(((-1)**(0.5) *x*p[1]))+ 1/(1/((-1)**(0.5) *x*p[2])))',
         '+ p[0]+ 1/((-1)**(0.5) *x*p[2])',
         '1/(+ 1/(p[0])+ 1/(1/((-1)**(0.5) *x*p[2])))',
         '+ p[0]+ ((-1)**(0.5) *x*p[1])',
         '1/(+ 1/(p[0])+ 1/(((-1)**(0.5) *x*p[1])))',
         '+ ((-1)**(0.5) *x*p[1])+ 1/((-1)**(0.5) *x*p[2])',
         '1/(+ 1/(((-1)**(0.5) *x*p[1]))+ 1/(1/((-1)**(0.5) *x*p[2])))']

j = (-1)**(0.5)

#setup variance for parameters
sig_v_o = 0.1
sig_v_f = 0.1
sig_r = 0.5
sig_f = 1
mega = {}
for fil in fils:
    
    script = open(fil,'r')
    long_string = script.read()
    alpha = long_string.split("\n")
    f = np.array([])
    v_o = np.array([])
    v_f = np.array([])
    phase = np.array([])
    for stg in alpha:
        g = stg.split(" ")
        w,x,y,z = float(g[0]),float(g[1]),float(g[2]),float(g[3])
        f = np.append(f,[w])
        v_o = np.append(v_o,[x])
        v_f = np.append(v_f,[y])
        phase = np.append(phase,[z])

    phase = phase / 180
    r = 2.99 * 1000
    I_f = v_f/(r)
    Z = v_o/I_f
    dx_1 = r/v_f
    dx_2 = v_o*r/(v_f**2)
    dx_3 = v_o/v_f
    sig_Z = (((dx_1*sig_v_o)**2) + ((dx_2*sig_v_f)**2) + ((dx_3*sig_r)**2))**0.5
    sig_th = np.full(len(f),0.9/180)
    smooth_f = np.linspace(np.min(f),np.max(f),1000)
    
    fntn_val = {}
    def fitfunction_a(x, test, *p):
            return np.absolute(eval(test))
        
    def fitfunction_b(x, test, *p):
            return np.angle(eval(test), deg=True)/180
    for test in fntns:
        if fil == 'boxf.txt':
            a, b, c = 2.32, 1e-7, -5.48e-05
            if bool('p[0]' in test) == False:
                a = 0
            if bool('p[1]' in test) == False:
                b = 0
            if bool('p[2]' in test) == False:
                c = 0.
        elif fil == 'boxf_1.text':
            a, b, c = 2.32, 1e-7, -5.48e-05
            if bool('p[0]' in test) == False:
                a = 0
            if bool('p[1]' in test) == False:
                b = 0
            if bool('p[2]' in test) == False:
                c = 0
        elif fil == 'boxa.text':
            a, b, c = 841, 1e-7, -3.95e-05
            if bool('p[0]' in test) == False:
                a = 0
            if bool('p[1]' in test) == False:
                b = 0
            if bool('p[2]' in test) == False:
                c = 0
        elif fil == 'boxc.text':
            a, b, c = 9.76, -7.5121, -4.505
            if bool('p[0]' in test) == False:
                a = 0
            if bool('p[1]' in test) == False:
                b = 0
            if bool('p[2]' in test) == False:
                c = 0
        guessparams = np.array([a, b, c])
        
        def fitfunction_1(x, *p):
            return fitfunction_a(x, test, *p)
        def fitfunction_2(x, *p):
            return fitfunction_b(x, test, *p)

        try:
            popt1, pcov1 = opt.curve_fit(fitfunction_1, f, Z, sigma=sig_Z, absolute_sigma=1, p0 = guessparams)
            for i in range(0,len(popt1)):
               print('param ',i,' = ',popt1[i], '+/-', np.sqrt(pcov1[i,i]))
            popt2, pcov2 = opt.curve_fit(fitfunction_2, f, phase, sigma=sig_th, absolute_sigma=1, p0 = guessparams)
            for i in range(0,len(popt2)):
               print('param ',i,' = ',popt2[i], '+/-', np.sqrt(pcov2[i,i]))
        except RuntimeError:
            print('fail')
            pass
        
        def bestfit_param(fitfunction, x, y, popt, dy):
            yfit= fitfunction(x,*popt)
            chisq = sum( (y-yfit)**2 / dy**2 )
            dof = len(x) - len(popt)
            return stat.chi2.cdf(chisq,dof)
        
        chrs = (bestfit_param(fitfunction_1, f, Z, popt1, sig_Z), bestfit_param(fitfunction_2, f, phase, popt2, sig_th))
        if chrs[0] < 1.:
            if chrs[1] < chrs[0]:
                fntn_val[test] = (popt2, chrs)
            else:
                fntn_val[test] = (popt1, chrs)
        elif chrs[1] <1 :
            fntn_val[test] = (popt2, chrs)
        
        print('confidence to reject this: ' + str((bestfit_param(fitfunction_1, f, Z, popt1, sig_Z), bestfit_param(fitfunction_2, f, phase, popt2, sig_th))))
        print('Done with the test of ' + str(test))
    print('done with the data of ' + str(fil))
    
    
    
    fig, ax = plt.subplots(2, sharex=True)
    ax[0].errorbar(f,Z, yerr= sig_Z, xerr = sig_f, fmt='.')
    ax[0].set_title('Impedance vs Freq of ' + str(fil))
    ax[1].errorbar(f,phase, yerr = sig_th, xerr=sig_f, fmt='.')
    ax[1].set_title('Arguement of Impedance/pi vs Freq of ' + str(fil))
    for fits in fntn_val.keys():
        ax[0].plot(smooth_f, fitfunction_a(smooth_f, fits, *fntn_val[fits][0]), '--')
        ax[1].plot(smooth_f, fitfunction_b(smooth_f, fits, *fntn_val[fits][0]), '--')
    plt.show()
    mega[fil] = str(fntn_val)
