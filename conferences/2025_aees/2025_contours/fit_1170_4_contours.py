import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse, Circle
from matplotlib.path import Path
import shapefile
from shapely.geometry import Point, Polygon
import numpy as np
from mapping_tools import get_field_data
mpl.style.use('classic')

# this is just for getting a test ellipse
def get_ellipse_path(lon, lat, ex, ey, angle):
    ellipse = Ellipse((lon, lat), width=ex, height=ey, angle=angle)
    
    # Create 20-point unit circle
    theta = np.linspace(0, 2 * np.pi, 20, endpoint=True)
    circle = np.column_stack([np.cos(theta), np.sin(theta)])

    # Get the path
    epath = ellipse.get_path()
    # Get the list of path vertices
    evertices = epath.vertices.copy()
    # Transform the vertices so that they have the correct coordinates
    evertices = ellipse.get_patch_transform().transform(circle)
    
    return evertices

evertices = get_ellipse_path(146.36, -38.41, 0.65, 0.4, 10.0)


def fitEllipse(cont,method):

    x=cont[:,0]
    y=cont[:,1]

    x=x[:,None]
    y=y[:,None]

    D=np.hstack([x*x,x*y,y*y,x,y,np.ones(x.shape)])
    S=np.dot(D.T,D)
    C=np.zeros([6,6])
    C[0,2]=C[2,0]=2
    C[1,1]=-1
    E,V=np.linalg.eig(np.dot(np.linalg.inv(S),C))

    if method==1:
        n=np.argmax(np.abs(E))
    else:
        n=np.argmax(E)
    a=V[:,n]

    #-------------------Fit ellipse-------------------
    b,c,d,f,g,a=a[1]/2., a[2], a[3]/2., a[4]/2., a[5], a[0]
    num=b*b-a*c
    cx=(c*d-b*f)/num
    cy=(a*f-b*d)/num

    angle=0.5*np.arctan(2*b/(a-c))*180/np.pi
    up = 2*(a*f*f+c*d*d+g*b*b-2*b*d*f-a*c*g)
    down1=(b*b-a*c)*( (c-a)*np.sqrt(1+4*b*b/((a-c)*(a-c)))-(c+a))
    down2=(b*b-a*c)*( (a-c)*np.sqrt(1+4*b*b/((a-c)*(a-c)))-(c+a))
    a=np.sqrt(abs(up/down1))
    b=np.sqrt(abs(up/down2))

    #---------------------Get path---------------------
    ell=Ellipse((cx,cy),a*2.,b*2.,angle)
    #ell_coord=ell.get_verts()
    
    # Create 100-point unit circle
    theta = np.linspace(0, 2 * np.pi, 100, endpoint=True)
    circle = np.column_stack([np.cos(theta), np.sin(theta)])

    # Get the path
    epath = ell.get_path()
    # Get the list of path vertices
    evertices = epath.vertices.copy()
    # Transform the vertices so that they have the correct coordinates
    evertices = ell.get_patch_transform().transform(circle)
       
    params=[cx,cy,a,b,angle]
    print(params)

    return params, evertices

def plotConts(contour_list):
    '''Plot a list of contours'''
    import matplotlib.pyplot as plt
    fig=plt.figure()
    ax2=fig.add_subplot(111)
    for ii,cii in enumerate(contour_list):
        x=cii[:,0]
        y=cii[:,1]
        ax2.plot(x,y,'-', label=str(ii))
    plt.legend()
    plt.show(block=False)

#-------------------Read in data-------------------

#params1,ell_fit=fitEllipse(evertices,1)
#params2,ell2=fitEllipse(evertices,2)

#plotConts([evertices,ell_fit])

################################################################################
# parse shapefile
################################################################################

inshp = '1170.4_2025_PGA-0.033_lev_nat_contours_trimmed.shp'

sf = shapefile.Reader(inshp)

shapes = sf.shapes()
levels = get_field_data(sf, 'LEVELS', 'float')

################################################################################
# fit contours
################################################################################

fitted_contours = []
for shape, level in zip(shapes, levels):
    # convert tuple array to numpy
    listarray = []
    for pp in shape.points:
        listarray.append([pp[0], pp[1]])
    points = np.array(listarray)
    
    method = 1
    params, ellipse_fit = fitEllipse(points, method)
    
    fitted_contours.append(ellipse_fit)


################################################################################
# write shapefile
################################################################################

# set shapefile to write to 
outshp = 'AS1170_4_fitted_contours.shp'
w = shapefile.Writer(outshp)
w.field('LEVELS','F', 5, 3)

# loop through fitted_contours
for fitted_contour, level in zip(fitted_contours, levels):
    new_contour = []
    for vert in fitted_contour:
        point = Point(vert)
        new_contour.append(vert)
        
    w.record(level)
    w.line([fitted_contour])

w.close()

# write projection file
prjfile = outshp.strip().split('.shp')[0]+'.prj'
f = open(prjfile, 'w')
f.write('GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378137.0,298.257223563]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]]')
f.close()
