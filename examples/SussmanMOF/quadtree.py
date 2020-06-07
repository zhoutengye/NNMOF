import numpy as np
from itertools import permutations, repeat
import matplotlib.pyplot as plt
import math

class Points():
    vertice = []
    centers = []
    vertice_distance = []
    centers_distance = []
    vertice_flag = []
    centers_flag = []
    center_account = 0
    vertex_account = 0

    def __init__(self,x_vertice,y_vertice):
        for i  in range(len(x_vertice)):
            for j in range(len(y_vertice)):
                self.add_vertex([x_vertice[i],y_vertice[j]])
        for i  in range(len(x_vertice)-1):
            for j in range(len(y_vertice)-1):
                self.add_center([(x_vertice[i]+x_vertice[i+1])/2,(y_vertice[j]+y_vertice[j+1])/2])

    def add_center(self,p):
        self.centers.append(p)
        self.centers_distance.append(distance(shape_type, shape_para,p).value)
        if (distance(shape_type, shape_para,p).value >= 0):
            self.centers_flag.append(1)
        else:
            self.centers_flag.append(0)
        self.center_account = self.center_account + 1

    def add_vertex(self,p):
        self.vertice.append(p)
        self.vertice_distance.append(distance(shape_type, shape_para,p).value)
        if (distance(shape_type, shape_para,p).value >=0):
            self.vertice_flag.append(1)
        else:
            self.vertice_flag.append(0)
        self.vertex_account = self.vertex_account + 1


class distance():

    def __init__(self,types,para,point):
        self.x    = point[0]
        self.y    = point[1]
        self.para = para
        if types == 'zalesak':
            self.value = self.zalesak()
        if types == 'circle':
            self.value = self.circle()

    def zalesak(self):
        x_center = self.para[0]
        y_center = self.para[1]
        radius   = self.para[2]
        notch_width  = self.para[3]
        notch_height = self.para[4]
        d1 = math.sqrt( ( self.x - self.para[0] )**2.0 + ( self.y - self.para[1] )**2.0 ) - self.para[2]
        d2 = - max( 0.5 -self.para[3]/2-self.x , self.x-0.5 - self.para[3]/2 )
        d3 = self.para[4] - self.y 
        return(- max(d1,min(d2,d3)))

    def circle(self):
        x_center = self.para[0]
        y_center = self.para[1]
        radius   = self.para[2]
        return( - math.sqrt( ( self.x - self.para[0] )**2.0 + ( self.y - self.para[1] )**2.0 ) + self.para[2])



class Cell():
    IN  =  1
    OUT = -1
    INTERFACE = 0
    LEAF = 1
    BRANCH = 0

    def __init__(self,center,rect,parent,level):
        self.center   = center
        self.parent   = parent
        self.level    = level
        self.children = [None,None,None,None]
        self.vertice  = rect
        self.get_status(rect)

    def get_status(self,rect):
        p = points.vertice_flag
        i1,i2,i3,i4 = rect
        distance_flag = [p[i1],p[i2],p[i3],p[i4]]
        if distance_flag == [1,1,1,1]:
            self.status = self.IN
            self.type   = self.LEAF
        elif distance_flag == [0,0,0,0]:
            self.status = self.OUT
            self.type   = self.LEAF
        else:
            self.status = self.INTERFACE
            self.type   = self.BRANCH

    def printout(self):
        print(['center   = ',self.center])
        if self.parent == None:
            print('parent   = None')
        else:
            print(['parent   = ',self.parent.center])
        print(['level    = ',self.level])
        if self.children[0] == None:
            print(['children = ',self.children])
        else:
            print(['children = ',self.children[0].center,self.children[1].center,self.children[2].center,self.children[3].center])
        print(['vertice  = ',self.vertice])
        print(['status   = ',self.status])
        print(['type     = ',self.type])
        print('')



class Tree():
    leaves = []
    all_cells = []
    level_size = []
    level_length = [0]
    rects = []
    max_depth = 0

    def __init__(self, x_vertice, y_vertice,maxdepth):
        n = 0
        depth = 0
        self.max_depth = maxdepth
        self.level_size.append([x_vertice[1]-x_vertice[0],y_vertice[1]-y_vertice[0]])
        for i  in range(len(x_vertice)-1):
            for j in range(len(y_vertice)-1):
                rect = []
                rect.append(points.vertice.index([x_vertice[i], y_vertice[j]])) 
                rect.append(points.vertice.index([x_vertice[i], y_vertice[j+1]])) 
                rect.append(points.vertice.index([x_vertice[i+1], y_vertice[j]])) 
                rect.append(points.vertice.index([x_vertice[i+1], y_vertice[j+1]])) 
                self.rects.append(rect)
                self.all_cells.append(Cell(n, rect, None, depth))
                n = n+1
        self.level_length.append((len(x_vertice)-1)*(len(y_vertice)-1))
        self.traverse(depth)
        self.find_leaves()
        self.quad_volume_fraction()


    def traverse(self, depth):
        depth = depth + 1
        self.level_length.append(self.level_length[depth])
        print(self.level_length, depth)

        if (depth >= self.max_depth+1):
            return
        else:
            self.level_size.append([self.level_size[depth-1][0]/2,self.level_size[depth-1][1]/2])
            for i in range(self.level_length[depth-1],self.level_length[depth]):   
                if self.all_cells[i].type == Cell.BRANCH:
                    self.subdevide(self.all_cells[i],depth)
                    self.level_length[depth+1] = self.level_length[depth+1] + 4
            self.traverse(depth)

    def subdevide(self,parent,depth):
        base_center  = points.centers[parent.center]
        base_vertex = [points.vertice[parent.vertice[0]],\
                        points.vertice[parent.vertice[1]],\
                        points.vertice[parent.vertice[2]],\
                        points.vertice[parent.vertice[3]]]
        sub_center = []
        sub_center.append([(base_center[0]+base_vertex[0][0])/2,(base_center[1]+base_vertex[0][1])/2])
        sub_center.append([(base_center[0]+base_vertex[1][0])/2,(base_center[1]+base_vertex[1][1])/2])
        sub_center.append([(base_center[0]+base_vertex[2][0])/2,(base_center[1]+base_vertex[2][1])/2])
        sub_center.append([(base_center[0]+base_vertex[3][0])/2,(base_center[1]+base_vertex[3][1])/2])
        sub_vertex = []
        sub_vertex.append([(base_vertex[0][0]+base_vertex[1][0])/2,(base_vertex[0][1]+base_vertex[1][1])/2])
        sub_vertex.append([(base_vertex[0][0]+base_vertex[2][0])/2,(base_vertex[0][1]+base_vertex[2][1])/2])
        sub_vertex.append([(base_vertex[1][0]+base_vertex[3][0])/2,(base_vertex[1][1]+base_vertex[3][1])/2])
        sub_vertex.append([(base_vertex[2][0]+base_vertex[3][0])/2,(base_vertex[2][1]+base_vertex[3][1])/2])
        sub_vertex.append(base_center)
        for i in range(4):
            points.add_center(sub_center[i])
        for j in range(5):
            if (sub_vertex[j] not in points.vertice):
                points.add_vertex(sub_vertex[j])

        rect = []
        rect.append(points.vertice.index(base_vertex[0])) 
        rect.append(points.vertice.index(sub_vertex[0])) 
        rect.append(points.vertice.index(sub_vertex[1]))
        rect.append(points.vertice.index(sub_vertex[4])) 

        rect.append(points.vertice.index(sub_vertex[0])) 
        rect.append(points.vertice.index(base_vertex[1])) 
        rect.append(points.vertice.index(sub_vertex[4]))
        rect.append(points.vertice.index(sub_vertex[2])) 

        rect.append(points.vertice.index(sub_vertex[1])) 
        rect.append(points.vertice.index(sub_vertex[4])) 
        rect.append(points.vertice.index(base_vertex[2]))
        rect.append(points.vertice.index(sub_vertex[3])) 

        rect.append(points.vertice.index(sub_vertex[4])) 
        rect.append(points.vertice.index(sub_vertex[2])) 
        rect.append(points.vertice.index(sub_vertex[3]))
        rect.append(points.vertice.index(base_vertex[3])) 
        
        current_account = len(points.centers)
        for i in range(4):
            sub_rect = rect[(i+1)*4-4:(i+1)*4]
            self.rects.append(sub_rect)
            self.all_cells.append(Cell(current_account-4+i, sub_rect, parent, depth))
            parent.children[i] = self.all_cells[current_account-4+i]

    def find_leaves(self):
        depth = self.max_depth + 1
        for i in range(self.level_length[depth-1],self.level_length[depth]):
        	if (self.all_cells[i].type == Cell.BRANCH):
        		self.all_cells[i].type = Cell.LEAF
        		# self.all_cells[i].printout()
            # print(i)
        for i in range(len(self.all_cells)):
        	if (self.all_cells[i].type == Cell.LEAF):
        		self.leaves.append(self.all_cells[i])



    def printtree(self):
        print(['cell_account = ', len(self.all_cells)])
        print(['leaf_account = ', len(self.leaves)])
        print(['level_size   = ', self.level_size])
        print(['level_length = ', self.level_length])
        print(['max_depth    = ', self.max_depth])


    def quadplot(self,x_vertice,y_vertice,plot_depth):

        X=[]
        Y=[]
        depth = 0
        # for i in range(len(x_vertice)):
        #     X.append([x_vertice[i],x_vertice[i]])
        #     Y.append([y_vertice[0],y_vertice[len(y_vertice)-1]])
        # for i in range(len(y_vertice)):
        #     X.append([x_vertice[0],x_vertice[len(y_vertice)-1]])
        #     Y.append([y_vertice[i],y_vertice[i]])

        depth = 1
        for n in range(plot_depth):
            for i in range(self.level_length[depth-1],self.level_length[depth]):
                if self.all_cells[i].type == Cell.BRANCH:
                    nn = self.all_cells[i].center
                    xc = points.centers[nn][0]
                    yc = points.centers[nn][1]
                    x0 = xc - self.level_size[depth-1][0]/2
                    x1 = xc + self.level_size[depth-1][0]/2
                    y0 = yc - self.level_size[depth-1][1]/2
                    y1 = yc + self.level_size[depth-1][1]/2
                    X.append([x0,x1])
                    X.append([xc,xc])
                    Y.append([yc,yc])
                    Y.append([y0,y1])
            depth = depth+1

        for i in range(len(X)):
            plt.plot(X[i], Y[i],color='k')

        plt.show()

    def quadplot2(self,x_vertice,y_vertice,plot_depth):

        X=[]
        Y=[]
        depth = 0
        for i in range(len(x_vertice)):
            X.append([x_vertice[i],x_vertice[i]])
            Y.append([y_vertice[0],y_vertice[len(y_vertice)-1]])
        for i in range(len(y_vertice)):
            X.append([x_vertice[0],x_vertice[len(y_vertice)-1]])
            Y.append([y_vertice[i],y_vertice[i]])

        depth = 1
        for n in range(plot_depth):
            for i in range(self.level_length[depth-1],self.level_length[depth]):
                if self.all_cells[i].type == Cell.BRANCH:
                    nn = self.all_cells[i].center
                    xc = points.centers[nn][0]
                    yc = points.centers[nn][1]
                    x0 = xc - self.level_size[depth-1][0]/2
                    x1 = xc + self.level_size[depth-1][0]/2
                    y0 = yc - self.level_size[depth-1][1]/2
                    y1 = yc + self.level_size[depth-1][1]/2
                    X.append([x0,x1])
                    X.append([xc,xc])
                    Y.append([yc,yc])
                    Y.append([y0,y1])
            depth = depth+1

        for i in range(len(X)):
            plt.plot(X[i], Y[i],color='k')

        plt.show()

    def quad_volume_fraction(self):
        self.volume = [0]*self.level_length[1]
        for i in range(len(self.leaves)):
            if self.leaves[i].status == 1 or self.leaves[i].status == 0:
                if self.leaves[i].level == 0:
                    n = self.leaves[i].center
                    self.volume[n] = 1
                else:
                    lvl = self.leaves[i].level
                    node = self.leaves[i].parent
                    n = node.center
                    lvl = lvl - 1   
                    while lvl>0:
                        node = self.all_cells[n].parent
                        n = node.center
                        lvl = lvl -1
                    self.volume[n] = self.volume[n] + pow(0.25,self.leaves[i].level)

    def quad_volume_plot(self,x_vertice,y_vertice):
        xc = np.zeros(len(x_vertice)-1)
        yc = np.zeros(len(x_vertice)-1)
        f  = np.zeros((len(x_vertice)-1,len(x_vertice)-1))
        nn = 0
        for i in range(len(x_vertice)-1):
            for j in range(len(y_vertice)-1):
                xc[i] = (x_vertice[i]+x_vertice[i+1])/2
                yc[j] = 0.5*(y_vertice[j]+y_vertice[j+1])
                f[i,j] = self.volume[nn]
                nn = nn+1
        Y, X = np.meshgrid(xc,yc)
        fig = plt.figure()
        plt.contour(X, Y, f,[0.5])
        plt.show()



    def write_quad_volume(self,x_vertice,y_vertice):
        f = open("DATA/quadtree.dat", "w")
        nn = 0
        for i in range(len(x_vertice)-1):
            for j in range(len(y_vertice)-1):
                xc = (x_vertice[i]+x_vertice[i+1])/2
                f.write('%12.8f'   %xc)
                yc = 0.5*(y_vertice[j]+y_vertice[j+1])
                f.write('%12.8f'   %yc)
                f.write('%12.8f\n' %self.volume[nn])
                nn = nn+1
        f.close()


global shape_type
global shape_para
global points

shape_type = 'zalesak'
shape_para = [0.5,0.75,0.15,0.1,0.85]

x_vertice = np.linspace(0,1,51)
y_vertice = np.linspace(0,1,51)
maxdepth = 3
points = Points(x_vertice,y_vertice)


b = Tree(x_vertice,y_vertice,maxdepth)
b.write_quad_volume(x_vertice,y_vertice)
b.quad_volume_plot(x_vertice,y_vertice)
b.quadplot(x_vertice,y_vertice,maxdepth)
# print(b.volume)