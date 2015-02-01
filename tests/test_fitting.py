if __name__ == '__main__':
    import os, commands
    from masai.api import Bruker, Process, Plot1d, Plot2d,  Fit

    # we try here a multi source fit
    # print(commands.getoutput('pwd'))

    def singlefit():

        DATADIR = os.path.join('..','doc','tutorial','exemples')
        user = 'user'
        name = 'SPEC1D'
        expno = 10

        path1 = os.path.join(DATADIR, user, 'nmr', name, str(expno))
        source = Bruker(path1)
        Process(source, transform="em 50; zf 1; ft --tdeff 512")
        Process(source, transform="pk --auto --fit_phc1 --bound_phc1=100; ab --mode poly")
        #p = Plot1d(source, xlim=(30, -20))

        # fit
        # try with variable width
        pos = [-94.9200018156, -95.0999986491, -94.3173715978, -93.9072128893,
               -92.8999998621, -92.6363220634, -82.0826836442]

        script = """

        COMMON:
        # common parameters ex.
        # $ gwidth: 1.0, 0.0, none
          $ gratio: 0.5, 0.0, 1.0

        MODEL: LINE_1
        shape: voigtmodel
            * ampl:  1.0, 0.0, none
            $ pos:   4.5, -100.0, 100.0
            > ratio: gratio
            $ width: 4.0, 0, 20

        MODEL: LINE_2
        shape: voigtmodel
            $ ampl:  0.2, 0.0, none
            $ pos:   1.6, 0, 2.0
            > ratio: gratio
            $ width: 4.0, 0, 20

        MODEL: LINE_3
        shape: voigtmodel
            $ ampl:  0.18, 0.0, none
            $ pos:   7.9, 5.0, 10.0
            > ratio: gratio
            $ width: 4.0, 0, 20

        MODEL: LINE_4
        shape: voigtmodel
            $ ampl:  0.1, 0.0, none
            $ pos:   14, 10.0, 20.0
            > ratio: gratio
            $ width: 6.0, 0, 20

        MODEL: LINE_5
        shape: voigtmodel
            $ ampl:  0.1, 0.0, none
            $ pos:   -5, -8.0, -1.0
            > ratio: gratio
            $ width: 6.0, 0, 10
        """
        f = Fit(source, script)

        f.run(maxiter=100, every=5, method='simplex')

        # fitted model plot
        p = Plot1d(source,
                    figsize=(3.4, 2.0), dpi=180,
                    lw=1, c='red',
                    modls = '--',
                    plot_model=True,
                    offset=15)

    def multifit():

        DATADIR = os.path.expanduser(os.path.join('..', '..', '..', 'exemples'))
        user = 'user'
        name = 'RUB-17'

        path1 = os.path.join(DATADIR, user, 'nmr', name, str(1))
        source1 = Bruker(path1)
        Process(source1, transform='em 1; zf 0; ft; cl -105 -75; pk --fit_phc1 --auto; ab --mode poly')

        path2 = os.path.join(DATADIR, user, 'nmr', name, str(2))
        source2 = Bruker(path2)
        Process(source2, transform='em 1; zf 0; ft; cl -105 -75; pk --fit_phc1 --auto; ab --mode poly') #')#

        #p = Plot1d(source1, hold=True)
        #p = Plot1d(source2, axe=p.axe)

        # fit
        # try with variable width
        pos = [-94.9200018156, -95.0999986491, -94.3173715978, -93.9072128893,
               -92.8999998621, -92.6363220634, -82.0826836442]

        script = """
        COMMON:
          experiment_variables: ampl
          $ gratio: 0.5, 0.0, 1.0
          $ gwidth: 0.6, 0.5, 0.7

        """
        for i, delta in enumerate(pos):
            fa = '*' if i==0 else '$'
            script +="""
                MODEL: LINE_%d
                shape: voigtmodel
                    %s ampl:  0.6, 0.5, 2
                    $ pos: %.3f, %.3f, %.3f
                    > ratio: gratio
                    > width: gwidth

            """%(i, fa, delta, delta-.5, delta+0.5)

        print(script)

        # prepare fit
        f = Fit( [source1, source2], script)

        # initial model plot
        #p = Plot1d([source1, source2],
        #            figsize=(3.4, 2.0), dpi=180,
        #            row = 2,
        #            col = 1,
        #            lw=.5,
        #            modls = '--',
        #            plot_model=True,
        #            offset=15,
        #            xlim=(-98,-80))

        f.run(maxiter=200, every=5, method='simplex')

        # fitted model plot
        p = Plot1d(source1,
                    figsize=(3.4, 2.0), dpi=180,
                    #ow = 2,
                    #col = 1,
                    lw=.5,
                    modls = '--',
                    plot_model=True,
                    offset=15,
                    xlim=(-98,-80))

            # fitted model plot
        p = Plot1d(source2,
                    figsize=(3.4, 2.0), dpi=180,
                    #ow = 2,
                    #col = 1,
                    lw=.5,
                    modls = '--',
                    plot_model=True,
                    offset=15,
                    xlim=(-98,-80))

    #### TEST 2D ######
    def spec2Q():

        # reading
        DATADIR = os.path.join( '..', '..', '..', 'exemples')
        user = 'user'
        name = 'RUB-17'
        expno = 102
        path = os.path.join(DATADIR, user, 'nmr', name, str(expno))
        source = Bruker(path)

        # scaling preparation
        # This is a DQ experiment so:
        source.par2.SFO1 = source.par.SFO1*2.0
        source.par2.SF = source.par.SF*2.0
        # correction due to the missetted ofset
        fn = 2.0 # 2 rotation of the spectra
        source.par2.SF = source.par2.SF + fn * source.par2.SW_h*1.0e-6

        # process  1srt dimention
        Process(source, "em 1")
        Process(source, "zf 2")
        Process(source, "ft")
        Process(source, transform = 'cl -99 -80')
        Process(source, 'pk 0 0 --auto --fit_phc1')
        Process(source, transform = 'ab --mode poly')

        # process 2nd dimension
        Process(source, axis=0, transform='em 1')
        Process(source, axis=0, transform='zf 2')
        Process(source, axis=0, transform='ft')
        Process(source, axis=0, transform = 'cl -198 -160')
        Process(source, axis=0, transform='pk 180 0 --auto --fit_phc1 --select cols')
        Process(source, axis=0, transform='ab --mode linear')

        # magnitude calculation
        Process(source, transform='mc')

        # plot2D
        p_2 = Plot2d(source, figsize=(3,2.7), spacing=5, exponent=1.2, xproj='sky', hold=True)

        # plot the diagonal
        p_2.axe.plot(source.data.columns, source.data.columns*2.0, lw=.1, ls='--')

        #row = source.row(-176.6, width=.2)
        #p_1 = Plot1d(row)

        common = "COMMON:"
        common +="\n * gwidth: 0.25, 0.21, 0.3"
        common +="\n * gwidth1: 1.5, 1.0, 2.0"

        spc={
            'A': -82.,
            'B': -92.6,
            'C': -92.7,
            'D': -93.9,
            'E': -94.2,
            'F': -94.9,
            'G': -95.0,
            }

        # make a dictionary sorted by key
        from collections import OrderedDict
        species = OrderedDict(sorted(spc.items(), key=lambda t: t[0]))
        print(species)
        connectivity = np.ones((7,7))  #*source.data.values.max()

        # some correlation are known to be missing
        connectivity[0,0] = 0.0                     # A-A missing
        connectivity[0,2] = connectivity[2,0] = 0.0  # A-C
        connectivity[0,3] = connectivity[3,0] = 0.0  # A-D
        connectivity[0,4] = connectivity[4,0] = 0.0  # A-E
        connectivity[0,5] = connectivity[5,0] = 0.0  # A-F

        connectivity = connectivity*source.data.values.max()

        def make_script(common, species, connectivity):
            script = common
            for i, key in enumerate(species.keys()):
                # additional common parameters
                v = species[key]
                dv = abs(v)*.01
                script += "\n $ gpos%s: %.3f, %.3f, %.3f"%(key, v, v-dv, v+dv)
            fixdone = False
            for i, key in enumerate(species.keys()):
                for j, key1 in enumerate(species.keys()):
                    # as it is normally symmetric we can try to reduce the number of varaible
                    if j<i:
                        continue
                    script += "\n\n MODEL: LINE_"+key+key1
                    script += "\n shape: gaussian2dmodel"
                    #if i==0:
                    #    script += "\n * ampl:  1., 0.0, none"
                    #else:
                    fi = '$' if connectivity[i,j]>0.001 else '*'
                    if fi=='$' and not fixdone:
                        fi = '*'
                        fixdone=True
                    script += "\n %s ampl:  %.3f, 0.0, none"%(fi, connectivity[i,j])
                    script += "\n > pos: gpos%s "%key
                    script += "\n > pos1: (gpos%s) + (gpos%s)"%(key, key1)
                    script += "\n > width: gwidth"
                    script += "\n > width1: gwidth1"
                    script += "\n * theta: 0.0"
            return script.lower()


        script = make_script(common, species, connectivity)
        print(script)
        f = Fit(source,  script, kind='SQ-DQ')

        f.run(maxiter=20000, maxfun=100000, every=100,
              constraints={'max_connections':3})
        x = source.data.columns
        y = source.data.index
        p_2.axe.contour(x, y, source.modeldata, 10)

        p_2.show()

        print(f.script)

    ###############################

    singlefit()
    #multifit()
    #spec2Q()