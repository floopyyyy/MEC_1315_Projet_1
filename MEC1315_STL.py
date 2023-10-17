"""
Couroisie de A. Comer sur StackOverflow pour le partage du code public
https://stackoverflow.com/questions/54006078/selection-of-face-of-a-stl-by-face-normal-value-threshold"
"https://stackoverflow.com/questions/59391743/how-to-generate-mesh-and-plot-3d-surface-in-python"
Remodifié par R.Nguyen pour adaptation du cours MEC1315
"""
import numpy as np
import struct
#####




def Unique(inputList):
      """ 
      Given an M x N list, this function gets the unique rows by treating all
      M Ntuples as single objects. This function also returns the indexing
      to convert the unique returned list back to the original non-unique list.
      """

      hashTable=dict()

      indexList=[]
      uniqueList=[]

      indx=0
      for ntuple in inputList:
            if not ntuple in hashTable:
                hashTable[ntuple]=indx
                indexList.append(indx)
                uniqueList.append(ntuple)
                indx+=1
            else:
                indexList.append(hashTable.get(ntuple))      

      return uniqueList, indexList


def IsBinarySTL(filename):
    try:
        with open(filename,'r') as f:
              test=f.readline()
    except UnicodeDecodeError:
        return True

    if len(test) < 5:
        return True
    elif test[0:5].lower() == 'solid':
        return False  # ASCII STL
    else:
        return True

def ReadSTL(filename):
    """ Returns numpy arrays for vertices and facet indexing """
    def GetListFromASCII(filename):
        """ Returns vertex listing from ASCII STL file """
        outputList=[]

        with open(filename,'r') as f:
            lines=[line.split() for line in f.readlines()]
        for line in lines:
            if line[0] == 'vertex':
                    outputList.append(tuple([float(x) for x in line[1:]]))
        return outputList

    def GetListFromBinary(filename):
        """ Returns vertex listing from binary STL file """
        outputList=[]
        with open(filename,'rb') as f:
            f.seek(80) # skip header
            nFacets=struct.unpack('I',f.read(4))[0] # number of facets in piece

            for i in range(nFacets):
                  f.seek(12,1) # skip normal
                  outputList.append(struct.unpack('fff',f.read(12))) # append each vertex triple to list (each facet has 3 vertices)
                  outputList.append(struct.unpack('fff',f.read(12))) 
                  outputList.append(struct.unpack('fff',f.read(12)))
                  f.seek(2,1) # skip attribute
        return outputList

    if IsBinarySTL(filename):
        vertexList = GetListFromBinary(filename)
    else:
        vertexList = GetListFromASCII(filename)

    coords, tempindxs = Unique(vertexList)

    indxs = list()
    templist = list()
    for i in range(len(tempindxs)):
        if (i > 0 ) and not (i % 3):
            indxs.append(templist)
            templist = list()
        templist.append(tempindxs[i])
    indxs.append(templist)

    return np.array(coords), np.array(indxs)


def GetNormals(vertices, facets):
    """ Returns normals for each facet of mesh """
    u = vertices[facets[:,1],:] - vertices[facets[:,0],:]
    v = vertices[facets[:,2],:] - vertices[facets[:,0],:]
    normals = np.cross(u,v)
    norms = np.sqrt(np.sum(normals*normals, axis=1))
    return normals/norms[:, np.newaxis]



def LireSTL(nom_fichier):
    vertex_fvn, face_fvn =ReadSTL(nom_fichier)
    normal_fvn=GetNormals(vertex_fvn, face_fvn)
    return face_fvn, vertex_fvn, normal_fvn


def Rx(angle):
    Rx=np.array([[1,     0      ,          0     ],
                 [0, np.cos(angle), np.sin(angle)],
                 [0, -np.sin(angle), np.cos(angle)]])
    return Rx

def Ry(angle):
    Ry=np.array([[np.cos(angle), 0, -np.sin(angle)],
                [0            , 1 ,     0    ],
                [np.sin(angle) ,0 ,    np.cos(angle)] ]  )
    return Ry

def Rz(angle):
    Rz=np.array([[np.cos(angle),  np.sin(angle),  0],
                 [-np.sin(angle), np.cos(angle),  0],
                 [0             , 0            ,  1]])
    return Rz


def EcrireSTLASCII(output_file_name,face_fvn,vertex_fvn,normal_fvn):
    face_fvn=face_fvn.astype(int)
    f=open(output_file_name,"w")
    f.write("solid python")
    for j in range(len(face_fvn)):
        t1="facet normal"
        t2="outer loop"
        t3="endloop"
        t4="endfacet"
        f.write("\n  %s  %7.6e  %7.6e  %7.6e" % (t1, normal_fvn[j,0],normal_fvn[j,1],normal_fvn[j,2]))    #Écrire facet normal
        f.write("\n    %s" % t2)   #Écrire outerloop
        f.write("\n      vertex  %7.6e %7.6e %7.6e" % (vertex_fvn[face_fvn[j,0],0], vertex_fvn[face_fvn[j,0],1], vertex_fvn[face_fvn[j,0],2])) #Écrire vertex
        f.write("\n      vertex  %7.6e %7.6e %7.6e" % (vertex_fvn[face_fvn[j,1],0], vertex_fvn[face_fvn[j,1],1], vertex_fvn[face_fvn[j,1],2])) #Écrire vertex
        f.write("\n      vertex  %7.6e %7.6e %7.6e" % (vertex_fvn[face_fvn[j,2],0], vertex_fvn[face_fvn[j,2],1], vertex_fvn[face_fvn[j,2],2])) #Écrire vertex
        f.write("\n    %s" % t3) #ecrire endloop
        f.write("\n  %s" % t4) #ecrire endfacet
    f.write("\n%s" % "Endsolid Python")
    f.close()    
    
def CalculNormal(face_fvn,vertex_fvn):
    face_fvn=face_fvn.astype(int)
    normal_fvn=np.empty(np.shape(face_fvn))
    for i in range(len(face_fvn)):
        # print(f[i])
        vecteurA=vertex_fvn[face_fvn[i],0]
        vecteurB=vertex_fvn[face_fvn[i],1]
        vecteurC=vertex_fvn[face_fvn[i],2]
        vecteur_normal=np.cross((vecteurB-vecteurA), (vecteurC-vecteurA))
        normal_fvn[i]=vecteur_normal/np.linalg.norm(vecteur_normal)
    return normal_fvn                


