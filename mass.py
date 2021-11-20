#!/usr/bin/env python
# coding: utf-8

# In[5]:


from pyopenms import *
seq = AASequence.fromString("VAKA")
V_weight=seq.getMonoWeight()
A_weight=seq.getMonoWeight()
K_weight=seq.getMonoWeight()
A_weight=seq.getMonoWeight()

print("Monoisotopic mass of peptide [V] is ",V_weight)
print("Monoisotopic mass of peptide [A] is ",A_weight)
print("Monoisotopic mass of peptide [K] is ",K_weight)
print("Monoisotopic mass of peptide [A] is ",A_weight)
print ("The piptide", str(seq), "consists of the following amino acids:")
for aa in seq:
    print(aa.getName(), ":", aa.getMonoWeight())


# In[ ]:





# In[ ]:




