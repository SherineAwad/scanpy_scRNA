with open(config['SAMPLES']) as fp:
    samples = fp.read().splitlines()
SUBSET = ['Cones', 'AC'] 

rule all:
         input:
            expand("{all}.h5ad", all= config['ALL']), 
            expand("doubletRemoved_{all}.h5ad", all=config['ALL']),
            expand("corrected_doubletRemoved_{all}.h5ad", all=config['ALL']),
            expand("clustered_corrected_doubletRemoved_{all}.h5ad", all=config['ALL']), 
            #expand("annotated_{all}.h5ad", all=config['ALL']), 
            #expand("{subset}"_{all}.h5ad", all=config['ALL'], subset = SUBSET),
 
rule preprocess: 
        input:  
            expand("{sample}_filtered_feature_bc_matrix.h5", sample = samples) 
        output: 
          expand("{all}.h5ad", all= config['ALL']), 
        params: 
          samplesFile = config['SAMPLES'],  
          name = config['ALL']
        shell: 
            """
           python preprocess.py {params.samplesFile}  {params.name}  
           """ 

rule remove_doublet: 
      input:
           expand("{all}.h5ad", all= config['ALL']),
      output:
          expand("doubletRemoved_{all}.h5ad", all= config['ALL']),
      params: 
         name = config['ALL'] 
      shell:
           """
           python doublet.py {input} {params.name}  
           """

rule batch: 
     input: 
         expand("doubletRemoved_{all}.h5ad", all=config['ALL'])
     output: 
         expand("corrected_doubletRemoved_{all}.h5ad", all=config['ALL'])
     params:
         name = config['ALL']
     shell: 
        """ 
        python batch.py {input} {params.name} 
        """ 


rule cluster: 
       input:
          expand("corrected_doubletRemoved_{all}.h5ad", all=config['ALL']) 
       output:
          expand("clustered_corrected_doubletRemoved_{all}.h5ad", all=config['ALL'])
       params:
         name = config['ALL']
       shell:
          """
          python cluster.py {input} {params.name} 
          """

rule annotate:
       input:
          expand("clustered_{all}.h5ad", all=config['ALL'])
       params: 
          annofile = config['ANNOFILE'], 
          name = config['ALL'] 
       output:
          expand("annotated_{all}.h5ad", all=config['ALL'])
       shell:
          """
          python annotate.py {input} {params.annofile} {params.name} 
          """

rule subset: 
      input: 
         expand("annotated_{all}.h5ad", all=config['ALL'])
      output: 
         expand("{subset}_{all}.h5ad", all=config['ALL'], subset = SUBSET),
      shell: 
          """ 
          python subset.py {input} ConesSubtypes.txt" 
          """
      
