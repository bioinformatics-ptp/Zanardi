#program for convert 705 interbull file to plink format
import os

def i705_convert_to_plink (gfile,mappa,outfile):    
    recode={'0':'B B','1':'A B','2':'A A','5':'0 0','7':'B B','8':'A B','9':'A A'}
    sexconv={'M':'1','F':'2','U':'0'}
    chrsex={'31':'X','32':'0','33':'Y','34':'MT'}
    dic_map={}

    for mapz in mappa:
        full = [line.strip().replace(";"," ").replace(',',' ').replace('\t',' ').split() for line in open(mapz)]
        dic_map[str(len(full))]=full
    
    for delfile in dic_map.keys():
        os.system('rm -f '+outfile+'_%s.ped' % (delfile))
        os.system('rm -f '+outfile+'_%s.map' % (delfile))

    error=False
    map_found={}
    for genofile in gfile:
        for line in open(genofile):
            data=line.strip().replace(";"," ").replace(',',' ').replace('\t',' ').split()
            info,idx,col,geno=data[:4]
            if dic_map.has_key(col) and int(col)==len(geno):
                if not map_found.has_key(col):map_found[col]=0
                map_found[col]+=1
                snp_recode=[recode[geno[x]] for x in range(len(geno))]
                o_pedfile=open(outfile+'_%s.ped' % col, 'a')
                try:o_pedfile.write('%s %s 0 0 %s -9 %s \n' % (idx[0:3],idx,sexconv[idx[6]],' '.join(snp_recode)))
                except KeyError: 
                    why_error={}
                    error=True
                    why_error[idx]=[col,idx[6],'sex']
                    break
                o_pedfile.close()                
            else:
                if not error: 
                    error = True
                    id_error=[]; why_error={}
                if int(col) != len(geno): why_error[idx]=[col,str(len(geno)),False]
                else: why_error[idx]=[col,str(len(geno)),True]

    map_warn={}
    for chips in map_found.keys():
        map_warn[chips]=False
        mp=open(outfile+'_%s.map' % chips, 'w')
        allSNP=dic_map[chips]
        for SNPinfo in allSNP:
            chromosome=SNPinfo[2]
            if SNPinfo[2] in chrsex.keys():
                chromosome=chrsex[SNPinfo[2]]
            ## ---- TODO GABRIELE: CONTROLLO NUMERO COLONNE (SE DIVERSO FERMATI)
            ## ---- FISSSA IL NUMERO DI COLONNE DI QUA SOTTO.
            if len(SNPinfo)==5: mp.write('%s %s %s %s \n' % ( chromosome,SNPinfo[0],'0',SNPinfo[3]))
            else:
                mp.write('%s %s %s %s \n' % ( SNPinfo[-3],SNPinfo[0],'0',SNPinfo[-1]))
                map_warn[chips]=True

    if error:return (False,map_found.keys(),map_found, map_warn, why_error)
    else:    return (True,map_found.keys(),map_found, map_warn)



