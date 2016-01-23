import subprocess as sub

def levels(structure, t):
  cmd = ['RNAshapes','-D', structure,'-t', str(t)]
  #cmd = ['RNAshapes','-D', structure,'-t', '%s' %t]
  p=sub.Popen(cmd,stdin=sub.PIPE,stdout=sub.PIPE,stderr=sub.PIPE)
  shapes_structure,error = p.communicate()
  shapes_structure = shapes_structure.replace('\n','')
  return shapes_structure
