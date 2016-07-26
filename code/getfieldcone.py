""" A script from Steve Myers """

# functions to pull out the fields from an MS that are within
# some distance, and match target field names

def getfieldict(msfile=None):
    #
    try:
        tb.open(msfile+'/FIELD')
    except:
        print 'ERROR: could not open '+msfile+'/FIELD'
        return
    field_dirs=tb.getcol('PHASE_DIR')
    field_names=tb.getcol('NAME')
    tb.close()
    #
    (nd,ni,nf) = field_dirs.shape
    print 'Found '+str(nf)+' fields'
    #
    # compile field dictionaries
    ddirs={}
    flookup={}
    for i in range(nf):
        fra = field_dirs[0,0,i]
        fdd = field_dirs[1,0,i]
        rapos = qa.quantity(fra,'rad')
        decpos = qa.quantity(fdd,'rad')
        ral = qa.angle(rapos,form=["tim"],prec=9)
        decl = qa.angle(decpos,prec=10)
        fdir = me.direction('J2000',ral[0],decl[0])
        ddirs[i]={}
        ddirs[i]['dir']=fdir
        fn=field_names[i]
        ddirs[i]['name']=fn
        if flookup.has_key(fn):
            flookup[fn].append(i)
        else:
            flookup[fn]=[i]
    print 'Cataloged '+str(nf)+' fields'
    #
    return ddirs

def getfieldcone(msfile=None,distance='0deg',center=''):
    # center = 'J2000 10:00:28.6 2.12.21'
    # e.g. mydir=me.direction('J2000','10:00:28.6','2.12.21')
    #
    fieldlist = []
    try:
        qdist = qa.toangle(distance)
        qdeg = qa.convert(qdist,'deg')
        maxdeg = qdeg['value']
    except:
        print 'ERROR: cannot parse distance ',distance
        return
    try:
        tb.open(msfile+'/FIELD')
    except:
        print 'ERROR: could not open '+msfile+'/FIELD'
        return
    field_dirs=tb.getcol('PHASE_DIR')
    field_names=tb.getcol('NAME')
    tb.close()
    #
    (nd,ni,nf) = field_dirs.shape
    print 'Found '+str(nf)+' fields'
    #
    # compile field dictionaries
    ddirs={}
    flookup={}
    for i in range(nf):
        fra = field_dirs[0,0,i]
        fdd = field_dirs[1,0,i]
        rapos = qa.quantity(fra,'rad')
        decpos = qa.quantity(fdd,'rad')
        ral = qa.angle(rapos,form=["tim"],prec=9)
        decl = qa.angle(decpos,prec=10)
        fdir = me.direction('J2000',ral[0],decl[0])
        ddirs[i]={}
        ddirs[i]['dir']=fdir
        fn=field_names[i]
        ddirs[i]['name']=fn
        if flookup.has_key(fn):
            flookup[fn].append(i)
        else:
            flookup[fn]=[i]
    print 'Cataloged '+str(nf)+' fields'
    #
    # turn the center into a direction
    if center=='':
        dcenter = ddirs[0]['dir']
    else:
        clis = center.split()
        if len(clis)==1:
            s = clis[0]
            if s.isdigit:
                # assume its a field index
                ic = int(s)
                if ic<0 or ic>(nf-1):
                    print 'ERROR: field index '+s+' out of range 0-'+str(nf-1)
                    return
                else:
                    dcenter=ddirs[ic]['dir']
                    fn = field_names[ic]
            else:
                # treat as field name
                if flookup.has_key(s):
                    # pick first field matching name
                    fn = s
                    ic = flookup[fn][0]
                else:
                    print 'ERROR: field name '+s+' unknown'
                    return
            print 'Using field index '+str(ic)+' name '+fn
            cdir = ddirs[ic]['dir']
        elif len(clis)==3:
            # assume its a direction EPO,RA,DEC
            try:
                cdir = me.direction(clis[0],clis[1],clis[2])
            except:
                print 'ERROR: not a direction ',clis
                return
            print 'Using center direction '+center
        elif len(clis)==2:
            # assume its a direction RA,DEC
            try:
                cdir = me.direction('J2000',clis[0],clis[1])
            except:
                print 'ERROR: not a direction ',clis
                return
            print 'Using center direction '+center
        else:
            print 'ERROR: unable to understand center '+center
            return
    #
    # Construct offset separations
    print 'Looking for fields with maximum separation '+str(maxdeg)+' deg'
    for i in range(nf):
        dd = ddirs[i]['dir']
        sep = me.separation(cdir,dd)
        sepdeg = qa.convert(sep,'deg')
        offs = sepdeg['value']
        ddirs[i]['offset']=offs
        if offs<=maxdeg:
            fieldlist.append(i)
    print 'Found '+str(len(fieldlist))+' fields within '+str(maxdeg)+' degrees'
    #
    return fieldlist

def getfieldircone(msfile=None,distance='0deg',center_dir=None):
    # Version STM 2016-May-16 use center direction measure
    #
    fieldlist = []
    try:
        qdist = qa.toangle(distance)
        qdeg = qa.convert(qdist,'deg')
        maxdeg = qdeg['value']
    except:
        print 'ERROR: cannot parse distance ',distance
        return
    try:
        tb.open(msfile+'/FIELD')
    except:
        print 'ERROR: could not open '+msfile+'/FIELD'
        return
    field_dirs=tb.getcol('PHASE_DIR')
    field_names=tb.getcol('NAME')
    tb.close()
    #
    (nd,ni,nf) = field_dirs.shape
    print 'Found '+str(nf)+' fields'
    #
    # compile field dictionaries
    ddirs={}
    flookup={}
    for i in range(nf):
        fra = field_dirs[0,0,i]
        fdd = field_dirs[1,0,i]
        rapos = qa.quantity(fra,'rad')
        decpos = qa.quantity(fdd,'rad')
        ral = qa.angle(rapos,form=["tim"],prec=9)
        decl = qa.angle(decpos,prec=10)
        fdir = me.direction('J2000',ral[0],decl[0])
        ddirs[i]={}
        ddirs[i]['dir']=fdir
        fn=field_names[i]
        ddirs[i]['name']=fn
        if flookup.has_key(fn):
            flookup[fn].append(i)
        else:
            flookup[fn]=[i]
    print 'Cataloged '+str(nf)+' fields'
    #
    # Construct offset separations
    print 'Looking for fields with maximum separation '+str(maxdeg)+' deg'
    for i in range(nf):
        dd = ddirs[i]['dir']
        sep = me.separation(center_dir,dd)
        sepdeg = qa.convert(sep,'deg')
        offs = sepdeg['value']
        ddirs[i]['offset']=offs
        if offs<=maxdeg:
            fieldlist.append(i)
    print 'Found '+str(len(fieldlist))+' fields within '+str(maxdeg)+' degrees'
    #
    return fieldlist

def getfieldirbox(msfile=None,distance='0deg',center_dir=None,matchregex=''):
    # Created STM 2016-May-16 use center direction measure
    # Returns list of fields from msfile within a 
    # rectangular box of size distance
    # Version STM 2016-Jun-07 add matchregex parameter for name
    import re
    #
    fieldlist = []
    #
    center_ra = center_dir['m0']['value']
    center_dec = center_dir['m1']['value']
    #
    try:
        qdist = qa.toangle(distance)
        qrad = qa.convert(qdist,'rad')
        maxrad = qrad['value']
    except:
        print 'ERROR: cannot parse distance ',distance
        return
    #
    try:
        tb.open(msfile+'/FIELD')
    except:
        print 'ERROR: could not open '+msfile+'/FIELD'
        return
    field_dirs=tb.getcol('PHASE_DIR')
    field_names=tb.getcol('NAME')
    tb.close()
    #
    (nd,ni,nf) = field_dirs.shape
    print 'Found '+str(nf)+' fields'
    #
    # compile field dictionaries
    ddirs={}
    flookup={}
    for i in range(nf):
        fra = field_dirs[0,0,i]
        fdd = field_dirs[1,0,i]
        rapos = qa.quantity(fra,'rad')
        decpos = qa.quantity(fdd,'rad')
        ral = qa.angle(rapos,form=["tim"],prec=9)
        decl = qa.angle(decpos,prec=10)
        fdir = me.direction('J2000',ral[0],decl[0])
        ddirs[i]={}
        ddirs[i]['ra']=fra
        ddirs[i]['dec']=fdd
        ddirs[i]['dir']=fdir
        fn=field_names[i]
        ddirs[i]['name']=fn
        if flookup.has_key(fn):
            flookup[fn].append(i)
        else:
            flookup[fn]=[i]
    print 'Cataloged '+str(nf)+' fields'
    #
    # Construct offset separations in ra,dec
    print 'Looking for fields with maximum separation '+distance
    nreject = 0
    skipmatch = matchregex=='' or matchregex==[]
    for i in range(nf):
        dd = ddirs[i]['dir']
        dd_ra = dd['m0']['value']
        dd_dec = dd['m1']['value']
        sep_ra_temp = abs(dd_ra - center_ra)
        if sep_ra_temp>(pl.pi):
            sep_ra_temp = 2.0*pl.pi - sep_ra_temp
        sep_ra = sep_ra_temp*pl.cos(center_dec)

        sep_dec = abs(dd_dec - center_dec)

        ddirs[i]['offset_ra']=sep_ra
        ddirs[i]['offset_ra']=sep_dec

        if sep_ra<=maxrad:
            if sep_dec<=maxrad:
                if skipmatch:
                    fieldlist.append(i)
                else:
                    # test regex against name
                    foundmatch = False
                    fn = ddirs[i]['name']
                    for rx in matchregex:
                        mat = re.findall(rx,fn)
                        if len(mat)>0:
                            foundmatch = True
                    if foundmatch:
                        fieldlist.append(i)
                    else:
                        nreject+=1
                        
    print 'Found '+str(len(fieldlist))+' fields within '+distance
    if not skipmatch:
        print 'Rejected '+str(nreject)+' distance matches for regex'
    #
    return fieldlist
