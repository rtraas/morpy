# This is so the package will monitor it's progress when running a Sim
#	Helps user debug
# Created by SF in March 2020

import time
import sys

class MonitorProgress:
    """
    Encapsulates the heartbeat needs for a rebound simulation
    only outputs if stdout or fileout is set to True

    Only outputs once per heartbeat operation and auto resets
    both the nextime and nextorbit if one triggers output

    Caviet Emptor
    No checks are done to make use time units are consitent between the
    simulation and internal time parameters
    """
       # timeinterval :   float       # number of seconds between outputs  
       # orbitperiod : float       # user defined length of an orbit
       # orbitinterval :  float       # number of orbits between outputs

       # stdout : bool             # output to standard output
       # fileout : bool            # output to a file
       # filename : str            # filename to output to

        # internalish
       # nextime = : float         # next time to output
       # nextorbit = : float       # next orbit to output
       # integrationtime : float         # the current time within the simulation
       # currentintegrationorbit : float # the current orbit in the simulation 
                                        # derived from user input orbit 
       # prestring : str           # text to prepend before output orbit
    

    def __init__(self,timeinterval=30, orbitinterval=10, orbitperiod=5000
                , stdout=False, fileout=True, filename='heartbeat', prestring=''):
        self.timeinterval=timeinterval         # number of seconds between outputs
        self.orbitperiod=orbitperiod           # length of an orbit
        self.orbitinterval=orbitinterval       # number of orbits between outputs

        self.stdout=stdout               # output to standard output
        self.fileout=fileout             # output to a file
        self.filename = filename         # filename to output to
        self.prestring = prestring       # text to prepend before output orbit 

        self.nextime = -1                # next time output
        self.nextorbit = -1              # next orbit output

    def heartbeat(self,rebsim):
        """ 
        pass to rebound.simuldation to monitor the progress of a long integration
        """
        self.currentorbit(rebsim)
        t=time.perf_counter()
        if t > self.nextime or self.currentintegrationorbit > self.nextorbit:
            self.updatenext(rebsim)     # update when to check next
            self.output(rebsim)         # output progress aka heartbeat 
        #print(f'{t:.2f} {self.nextime:.2f} {self.timeinterval} {self.orbitinterval} {self.currentintegrationorbit} {self.nextorbit}')


    def updatenext(self, rebsim):
        self.nextime= time.perf_counter()+self.timeinterval
        self.nextorbit= self.nextorbit+self.orbitinterval

    def output(self,rebsim):
        # define the output string
        outstring=f'{self.prestring}Current Orbit - {self.currentintegrationorbit:.2f}'
        # use rebsim.contents
        if self.stdout:
            print('\r'+outstring,end='',file=sys.__stdout__,flush=True)
        if self.fileout and self.filename != '':
            # open file to write - not append 
            # TODO move to a try block in case out directory is not accessible....
            f = open (self.filename,'w')
            # write file
            print(outstring,file=f,flush=True)
            # close file
            f.close()
    
    def currentorbit(self,rebsim):
        # figure out the orbit number of the current simuldation
        self.integrationtime= rebsim.contents.t
        self.currentintegrationorbit= self.integrationtime / self.orbitperiod
        return self.currentintegrationorbit

    def dictSet(self, inputdict):
        """
        set the properties using an input dictionary
        uses self reflection to only insert properities already defined in this class
        :param inputDict:  dictionary or properties to assign
        :return:
        """
        #    a=dir(self)
        #    print(a)

        # only keys beginning with hb_ are used

        #hbKeys =(k.lstrip('hb_') for k in inputdict.keys() if k.startswith('hb_'))

        relaventKeys = (k for k in inputdict.keys() \
            if k.startswith('hb_') and k[3:] in dir(self))

        for k in relaventKeys:
            #   print(k)
            setattr(self, k[3:], inputdict[k])   
    
    if __name__ == "__main__":
        import rebound
        import numpy as np
        import MonitorProgress

        # simple test of dictSet
        hb=MonitorProgress.MonitorProgress()
        testdict={
            "hb_timeinterval":4,
            "hb_filename":"goodfilename",
            "hb_orbitperiod":2.0,
            'filename':"badfilename",
            "orbitperiod":-0.1,
            "hb_orbitPeriod":-0.3,
            "timeinterval":"a",
            "hhbb__timeinterval":"b"
            }
        hb.dictSet(testdict)
        print(f'{hb.timeinterval} = 4 {hb.filename} = goodfilename {hb.orbitperiod} = 2.0')

        sim=rebound.Simulation()
        #print("hi")

        sim.dt=0.0001
        sim.integrator='ias15'

        #sim.add("Sun")
        #sim.add("Mercury")
        #sim.add("Venus")
        #sim.add("Earth")
        #sim.add("Mars")
        #sim.add("Jupiter")

        # add really slow satelites
        #sim.add("Callisto")
        #sim.add("Kale")
        #sim.add("Deimos")
        #sim.add("Phobos")

        #sim.units=('yr','AU','Msun')
        #sim.add(m=1.0, x=-0.0042115824023980395, y=0.007328697950525469, z=3.442380655558414e-05, vx=-0.0004787902669331758, vy=-0.00014296804757355929, vz=1.3269511890473582e-05)
        #sim.add(m=1.6601141530543488e-07, x=-0.18396228762368783, y=0.2750567575806898, z=0.038401023991617315, vx=-1.6882186703266158, vy=-0.8483962871912428, vz=0.08552301032494537)
        #sim.add(m=2.4478382877847715e-06, x=0.0713494892258671, y=0.7233442966392737, z=0.005498759377651357, vx=-1.1737838956321165, vy=0.11726674345071304, vz=0.06933225038688995)
        #sim.add(m=3.040432648022642e-06, x=-0.8632418549898224, y=0.49641724174338514, z=1.4198438013165557e-05, vx=-0.511652256270947, vy=-0.8730527069897175, vz=5.469367608945412e-05)
        #sim.add(m=3.2271560375549977e-07, x=-0.7998027219341177, y=-1.298402550273837, z=-0.007805913191787548, vx=0.7248628337448288, vy=-0.353630705150824, vz=-0.02518993830975503)
        #sim.add(m=0.0009547919152112404, x=0.890114057069815, y=-5.124688568728187, z=0.0013409032747576652, vx=0.4268282490040009, vy=0.09596835473764002, vz=-0.009946417664484843)
        #sim.add(m=5.40966985863744e-08, x=0.8995837617180917, y=-5.132880282112665, z=0.0012099202980875034, vx=0.6068056076310155, vy=0.30608472194307296, vz=-0.0008918841731685606)
        #sim.add(m=0.0, x=0.8315614702094507, y=-5.099210661194104, z=0.036192422377282416, vx=0.3838268001535819, vy=-0.008823162926506554, vz=0.01924650743472604)
        #sim.add(m=0.0, x=0.7177998760788925, y=-5.051974422836574, z=-0.003505001173678271, vx=0.4653362411401482, vy=0.14324213598818158, vz=0.005986079732203662)
        #sim.add(m=7.2454169670014585e-16, x=-0.7998519799952434, y=-1.2982579412555009, z=-0.007770641834123271, vx=0.6861122764265917, vy=-0.3707791909400132, vz=-0.008970484940618336)
        #sim.add(m=5.340528788902282e-15, x=-0.7998600263317945, y=-1.2984034989846664, z=-0.007778582449441869, vx=0.7238721425561384, vy=-0.424381411316606, vz=-0.028520026813855983)

        #sim.add(m=3.040432648022642e-28, x=-0.8632418749898224, y=0.49641724174338514, z=1.4198438013165557e-05, vx=-0.511652256270947, vy=-0.8730527069897175, vz=5.469367608945412e-05)
        #sim.add(m=3.040432648022642e-28, x=-0.8632418649898224, y=0.49641724174338514, z=1.4198438013165557e-05, vx=-0.511652256270947, vy=-0.8730527069897175, vz=5.469367608945412e-05)

        sim.units=('m','s','kg')
        sim.add(m=1.9884754159665356e+30, x=-630043759.6763372, y=1096357608.402064, z=5149728.162104089, vx=-2.2696911186239674, vy=-0.6777358067514885, vz=0.06290372918229457)
        sim.add(m=3.301096181046679e+23, x=-27520366517.60466, y=41147905255.71728, z=5744711421.845565, vx=-8002.950742668865, vy=-4021.7975408014427, vz=405.41930439786114)
        sim.add(m=4.867466257521636e+24, x=10673731663.722311, y=108210766560.2245, z=822602694.3883002, vx=-5564.287887815159, vy=555.899533682848, vz=328.6674851293504)
        sim.add(m=6.045825574495058e+24, x=-129139143405.59561, y=74262962343.57758, z=2124056.094035506, vx=-2425.4724084526197, vy=-4138.680570592591, vz=0.2592737560441174)
        sim.add(m=6.41712044416609e+23, x=-119648784181.4082, y=-194238256832.4157, z=-1167747992.360458, vx=3436.1908534804456, vy=-1676.3759127659516, vz=-119.41215853562089)
        sim.add(m=1.8985802507611564e+27, x=133159167617.7826, y=-766642497882.3674, z=200596274.7184038, vx=2023.3667074051236, vy=454.9351510676732, vz=-47.15069915645675)
        sim.add(m=1.0756995522395713e+23, x=134575815269.3227, y=-767867960762.07, z=181001500.3105998, vx=2876.544060081236, vy=1450.985583709955, vz=-4.227950579798423)
        sim.add(m=0.0, x=124399825299.49532, y=-762831057165.3771, z=5414309323.116482, vx=1819.5196092405722, vy=-41.82594324804678, vz=91.23750002049378)
        sim.add(m=0.0, x=107381333050.12619, y=-755764616487.2129, z=-524340712.38327026, vx=2205.9127067365043, vy=679.0351147854021, vz=28.376834162867077)
        sim.add(m=1440733351730922.0, x=-119656153082.4674, y=-194216623631.18863, z=-1162471472.357184, vx=3252.494981067886, vy=-1757.6678031439244, vz=-42.52431890457192)
        sim.add(m=1.0619510204993724e+16, x=-119657356797.2824, y=-194238398757.5357, z=-1163659371.500894, vx=3431.4945111619754, vy=-2011.767545618968, vz=-135.19834473026341)
        sim.add(m=604.5825574495058, x=-129149150905.55302, y=74262962343.57758, z=2124056.094035506, vx=0, vy=0, vz=0.2592737560441174)
        sim.add(m=604.5825574495058, x=-129150151001.5743, y=74262962343.57758, z=2124056.094035506, vx=0, vy=0, vz=0.2592737560441174)
        #sim.add(m=6.045825574495058e+24, x=-129,139,143,405.59561, y=74262962343.57758, z=2124056.094035506, vx=-2425.4724084526197, vy=-4138.680570592591, vz=0.2592737560441174)

        #sim.convert_particle_units('m','s','kg')

        p=sim.particles
        #p[11].x=p[11].x+5000
        #p[12].x=p[12].x+7000


        orbelements=p[11].calculate_orbit(primary=p[3])
        print(f'{orbelements.a}   {orbelements.P} {sim.G}')
        #orbelements=p[12].calculate_orbit(primary=p[3])
        #print(f' {orbelements.a}   {orbelements.P}')

        #sim.convert_particle_units=('yr','AU','Msun')

        #orbelements=p[11].calculate_orbit(primary=p[3])
        #print(f'{orbelements.a}   {orbelements.P}')

        sim.status()

        sim.move_to_com()

        #er=MonitorProgress()
        #mp=MonitorProgress.MonitorProgress(50000,6000,0.5,True)
        mp=MonitorProgress.MonitorProgress(5,7000,20000,True,True)
        sim.heartbeat=mp.heartbeat

        times = np.logspace(3,5,num=1000)
        for t in times:
            sim.integrate(t,1)
            #print(f'{sim.t}',sep=" ")


        

