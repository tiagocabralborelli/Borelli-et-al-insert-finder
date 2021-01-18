def create (status, path ):
    """-> Cria diretórios onde os resultados do ABRicate serão armazenados"""
    import os
    db = ['argannot', 'card', 'ncbi', 'plasmidfinder', 'resfinder', 'vfdb','ecoh','ecoli_vf']
    if status == 'yes':
        for databank in db:
            os.system(f'mkdir {path}/{databank}')
            
            

def execute (path, ext, status):
    import glob
    import os
    """-> Executa o ABRicate em loop para todos os arquivos num dado diretório"""
    banks = ['argannot', 'card', 'ncbi', 'plasmidfinder', 'resfinder', 'vfdb','ecoh','ecoli_vf']

    if status == 'yes':
        all_files = glob.glob(path + f"/*.{ext}")
        finalfiles = []
        if ext == 'gbff':
            print('abricate iniciado')
            for files in all_files:
                finalfiles.append(files[0:-5])   
            for files in finalfiles:
                for db in banks:
                    os.system(f'abricate --db {db} {files}.{ext} > {files}_{db}.csv')
                    os.system(f'mv {files}_{db}.csv {path}/{db}')
        elif ext == 'fasta':
            print('abricate iniciado')
            for files in all_files:
                finalfiles.append(files[0:-6])   
            for files in finalfiles:
                for db in banks:
                    os.system(f'abricate --threads 8 --db {db} {files}.{ext} > {files}_{db}.csv')
                    os.system(f'mv {files}_{db}.csv {path}/{db}')
        print('abricate finalizado')
        print(f'Foram analisados {len(all_files)} genomas')        
        
        
        
def concatcsv (path,database):
    """-> Agrupa todos os CSVs de um mesmo banco de dados resultantes do ABRIcate"""
    import pandas as pd
    import glob
    db =['argannot', 'card', 'ncbi', 'plasmidfinder', 'resfinder', 'vfdb','ecoh','ecoli_vf'] # list with all databases 
    
    #Concat
    path = f'{path}/{database}'

    all_files = glob.glob(path + "/*.csv")
    li = []
    for filename in all_files:
        df = pd.read_csv(filename, index_col=None, header=0, sep = '\t', names = ['path','sequence','start','end','strand','gene','coverage','coverage_map','gaps','coverage_perc','identity_perc','database','accession','product','resistance'])
        li.append(df)
    genomes = pd.concat(li, axis=0, ignore_index=True)
    genomes['file'] = genomes['path'].apply(lambda x: x.split('/')[-1])
    
    #Cuts
   # genomes.drop(genomes[genomes.gaps != '0/0'].index,axis = 0,inplace = True) #drops all lines with morer than 0 GAPS
    genomes.drop(genomes[genomes.coverage_perc < 80].index, axis = 0,inplace = True)# Drops all lines with percentual of coverage less than 100
   # genomes.drop(genomes[genomes.identity_perc < 50].index,axis = 0, inplace = True) # Drops all lines with percentual of identity less than 100
    
    #Reseting index
    genomes.reset_index(inplace = True)
    genomes.sort_values(by = ['start'])
    genomes.drop(['index'],axis = 1, inplace = True)
    return genomes

def getspecies(df):
    """-> Cria uma coluna com as espécies dos arquivos analisados. Apenas para aquivos GenBank"""
    from Bio import SeqIO
    import pandas as pd
    index = 0
    names_gen = []
    species = []
    
    for items in df['path'].iteritems():
        record = SeqIO.parse(f"{df.path[index]}", "genbank")
        for ind,info in enumerate(record):
            if ind == 0:
                names_gen.append(f"{info.annotations['source']} {df.file[index]}")
                index = index + 1

            else:
                break
    df['organism'] = names_gen
    df['specie'] =  df.organism.str.split(' ').apply(lambda x: x[0:2])
    df['specie'] = df.specie.apply(lambda x: " ".join(x[0:2]))
    return df

def specieshc(df, path):
    """-> Cria uma coluna com o local de coleta e especie dos genomas analisados. É especial para os genomas do hc"""
    import pandas as pd
    info = pd.read_csv(f'{path}/info.csv')
    info['assembly'] = info.assembly.apply(lambda x: f'{x}.fasta')
    isolation_source = []
    name_specie = []
    x = 0
    y = 0
    x,y
    
    while True:
        if df['file'][x] == info['assembly'][y]:
            isolation_source.append(info['source'][y])
            name_specie.append(info['specie'][y])
        y = y + 1
        if y == len(info):
            y = 0 
            x = x + 1
        if x == len(df):
            break
    df['source'] = isolation_source
    df['specie'] = name_specie
    return df
                                 
def getsource(df):
    """-> Adiciona uma coluna com a fonte de coleta da bactéria"""
    from Bio import SeqIO
    lista = []
    lista2 = []
    for items in df.path:
        for gb_record in SeqIO.parse(open(f"{items}","r"), "genbank"):
            for i,f in enumerate(gb_record.features):
                if i == 0:
                    try:
                        lista.append(f.qualifiers['isolation_source'][0])
                        file = open('/home/tiagoboreli/Documentos/abricate/brasil/errorfile.txt','a')
                        file.write(f'{items} \n')
                        file.close()
                    except:
                        lista.append('not available')
                        file = open('/home/tiagoboreli/Documentos/abricate/brasil/errorfile.txt','a')
                        file.write(f'{items} \n')
                        file.close()
        lista2.append(lista[0])
        lista = []
    df['source'] = lista2
    return df
                                 
def getlocation(df):
    """-> Adiciona uma coluna com o local de coleta """
    from Bio import SeqIO
    lista = []
    lista2 = []
    for items in df.path:
        for gb_record in SeqIO.parse(open(f"{items}","r"),"genbank"):
            for i,f in enumerate(gb_record.features):
                if i == 0:
                    try:
                        lista.append(f.qualifiers['country'][0])
                        file = open('/home/tiagoboreli/Documentos/abricate/brasil/errorfile.txt','a')
                        file.write(f'{df.path[index]} \n')
                        file.close()
                    except:
                        lista.append('not available')
        lista2.append(lista[0])
        lista = []
    df['location'] = lista2
    return df             
               
def pca (df):
    """-> Criar um dataframe para plotar um PCA"""
    from sklearn.preprocessing import StandardScaler
    from sklearn.decomposition import PCA
    import pandas as pd

    del df.index.name
    df.reset_index(inplace = True)
    df.rename(columns = {'index':'Grupos:'}, inplace = True)
    grupos = df['Grupos:']
    df.drop('Grupos:', inplace = True, axis = 1)
    genes = df.columns
    df['Grupos:'] = grupos
    features = list(genes)
    #Separar as caracteristicas
    x = df.loc[:, features].values
    #Separar os alvos
    y = df.loc[:,['Grupos:']].values
    #Normalização das caracteristicas
    x = StandardScaler().fit_transform(x)
    #PCA
    pca = PCA(n_components=2)
    principalComponents = pca.fit_transform(x)
    principalDf = pd.DataFrame(data = principalComponents, columns = ['PC1','PC2'])
    finalDf = pd.concat([principalDf, df[['Grupos:']]], axis = 1)
    finalDf.drop(0, axis = 0, inplace = True)                            
    
    
    return finalDf

def argextract(df):
    "Retira dos os genes encontrados para um dado dataset"
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.Alphabet import IUPAC
    df['gene'] = df['gene'].apply(lambda x: x.split(')')[1])
    df['gene'] = df['gene'].apply(lambda x: x.replace("''",''))
    df['gene'] = df['gene'].apply(lambda x: x.replace('(',''))
    
    count = 0
    for items in df.itertuples():
        for record in SeqIO.parse(f'{items.path}','fasta'):
            if record.id == items.sequence:
                if items.strand == '+':
                    i = items.start - 1  
                    f = items.end + 1
                    seq = record.seq
                    seq = seq[i:f]
                    dna = Seq(str(seq), IUPAC.unambiguous_dna)
                    count = count + 1
                    file = open('seqs.fasta','a')
                    file.write(f'>{items.specie}_{items.file[0]}{items.file[6:8]}_{items.gene}_{items.sequence}_{count}\n')
                    file.write(f'{dna}\n')
                    file.close()
                    
                elif items.strand == '-':
                    i = items.start
                    f = items.end 
                    seq2 = record.seq
                    seq2 = seq2[i:f]
                    dna2 = Seq(str(seq2), IUPAC.unambiguous_dna)
                    p_seq = dna2.reverse_complement()
                    count = count + 1
                    file = open('seqs.fasta','a')
                    file.write(f'>{items.specie}_{items.file[0]}{items.file[6:8]}_{items.gene}_{items.sequence}_{count}\n')
                    file.write(f'{p_seq}\n')
                    file.close()
    return print('file ready')
  
def betaextrac(df,gene,specie):
    """-> Encontra as sequencias dos ARGs encontrados """
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.Alphabet import IUPAC
    count = 0
    for items in df.itertuples():
        for record in SeqIO.parse(f'{items.path}','fasta'):
            if record.id == items.sequence:
                if items.strand == '+' and '(Bla)' in items.gene and items.specie == f'{specie}':
                    i = items.start - 1  
                    f = items.end + 1
                    seq = record.seq
                    seq = seq[i:f]
                    dna = Seq(str(seq), IUPAC.unambiguous_dna)
                    name = items.gene.split(')')[1]
                    count = count + 1
                    file = open(f'{gene}.fasta','a')
                    file.write(f'>{name}_{items.file[:-6]}_{count}\n')
                    file.write(f'{dna}\n')
                    file.close()
                    
                elif items.strand == '-' and '(Bla)' in items.gene and items.specie == f'{specie}':
                    i = items.start
                    f = items.end 
                    seq2 = record.seq
                    seq2 = seq2[i:f]
                    dna2 = Seq(str(seq2), IUPAC.unambiguous_dna)
                    p_seq = dna2.reverse_complement()
                    gene = items.gene.split(')')[1]
                    count = count + 1
                    file = open(f'{gene}.fasta','a')
                    file.write(f'>{gene}_{items.file[:-6]}_{count}\n')
                    file.write(f'{p_seq}\n')
                    file.close() 
    return print('file ready')