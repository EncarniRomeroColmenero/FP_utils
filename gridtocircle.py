#!/usr/bin/env python
import math
import numpy as np

#pixelarea returns the fractional area of a given pixel (pixX,pixY) 
#at a radius R from the central pixel (cirX,cirY)
def pixelarea(pixX,pixY,cirX,cirY,R):
    pixX=float(pixX)
    pixY=float(pixY)
    cirX=float(cirX)
    cirY=float(cirY)
    R=float(R)
#    R=round(R,3)
    # assumes pixX,pixY refer to center of pixel
    x1 = pixX - .5
    x2 = pixX + .5
    y1 = pixY - .5
    y2 = pixY + .5
    cornerUR = np.sqrt((x2-cirX)**2+(y2-cirY)**2)
    cornerLR = np.sqrt((x2-cirX)**2+(y1-cirY)**2)
    cornerUL = np.sqrt((x1-cirX)**2+(y2-cirY)**2)
    cornerLL = np.sqrt((x1-cirX)**2+(y1-cirY)**2)

    #comparisons simplified?
    lowerleft = (round(cornerLL,10)>round(R,10)) # (round(cornerLL,10)==round(R,10)) + 
    lowerright = (round(cornerLR,10)>round(R,10)) #  (round(cornerLR,10)==round(R,10)) +
    upperright = (round(cornerUR,10)>round(R,10)) # (round(cornerUR,10)==round(R,10)) +
    upperleft = (round(cornerUL,10)>round(R,10)) #  (round(cornerUL,10)==round(R,10)) +
    total = lowerleft + lowerright + upperright + upperleft
#    print "cornerLL, cornerLR, cornerUR, cornerUL"
#    print cornerLL, cornerLR, cornerUR, cornerUL
#    print "lowerleft, lowerright, upperright, upperleft, R, total"
#    print lowerleft, lowerright, upperright, upperleft, R, total
#    return ""
#    raise ValueError("dumb")
    
    if total >= 4:
        if x1 < cirX and x2 > cirX and y1 < cirY and y2 > cirY:
            return R**2 * math.pi
        else:
            return 0
    elif total == 3:
        if not(upperright):
            A=onecorner(x2,y2,cornerUR,cornerLR,cornerUL,cornerLL,R,cirX,cirY)
        elif not(upperleft):
            A=onecorner(x1,y2,cornerUR,cornerLR,cornerUL,cornerLL,R,cirX,cirY)
        elif not(lowerright):
            A=onecorner(x2,y1,cornerUR,cornerLR,cornerUL,cornerLL,R,cirX,cirY)
        elif not(lowerleft):
            A=onecorner(x1,y1,cornerUR,cornerLR,cornerUL,cornerLL,R,cirX,cirY)
        else:
            A=0
            raise ValueError("There was a problem with the corner - radius comparison")
    elif total == 2:
        A=twocorner(x1,x2,y1,y2,cornerUR,cornerLR,cornerUL,cornerLL,R,cirX,cirY)
    elif total == 1:
        if lowerleft:
            A=threecorner(x2,y2,x1,y1,R,cirX,cirY)
        elif lowerright:
            A=threecorner(x1,y2,x2,y1,R,cirX,cirY)
        elif upperleft:
            A=threecorner(x2,y1,x1,y2,R,cirX,cirY)
        elif upperright:
            A=threecorner(x1,y1,x2,y2,R,cirX,cirY)
        else:
            raise ValueError("Total=1 but other parameters not satisfied")
    elif total == 0:
        A=1
    else: 
        raise ValueError("I don't know what happened but total is %d" % total)
    if A > 1:
        raise ValueError("Area > 100%  %f" % A)
    elif A < 0:
        raise ValueError("Negative area?!")
    return A

def onecorner(x1,y1,cornerUR,cornerLR,cornerUL,cornerLL,R,cirX,cirY):
    if y1 > cirY: 
        y2 = cirY + np.sqrt( R**2 - (cirX-x1)**2 )
    elif y1 < cirY: 
        y2 = cirY - np.sqrt( R**2 - (cirX-x1)**2 )
    else:
        y2 = y1 + R**2
    if x1 > cirX:
        x2 = cirX + np.sqrt( R**2 - (cirY-y1)**2 )
    elif x1 < cirX:
        x2 = cirX - np.sqrt( R**2 - (cirY-y1)**2 )
    else:
        x2 = x1 + R**2
    if cornerUR < R:
        theta = abs( abs(math.acos( np.clip(abs(cirY - y2) / R,-1.0,1.0) )) - abs(math.asin ( np.clip(abs(cirX - x2) / R,-1.0,1.0) )) )
        t2 = math.pi/2 -  abs( abs(math.asin( np.clip(abs(cirY - y2) / R,-1.0,1.0) )) + abs(math.asin ( np.clip(abs(cirX - x2) / R,-1.0,1.0) )) )
#DEBUG  #print "Theta: %f 90-t2: %f" % (theta,t2)
    elif cornerLR < R:
        theta = abs( abs(math.asin( np.clip(abs(cirY - y1) / R,-1.0,1.0) )) - abs(math.asin ( np.clip(abs(cirY - y2) / R,-1.0,1.0) )) )
        t2 = math.pi/2 -  abs( abs(math.asin( np.clip(abs(cirY - y1) / R,-1.0,1.0) )) + abs(math.asin ( np.clip(abs(cirX - x2) / R,-1.0,1.0) )) )
#DEBUG  #print "Theta: %f 90-t2: %f" % (theta,t2)
#DEBUG  #print "Corner 2: %f Theta: %f x2: %f y2: %f"% (cornerLR,theta,x2,y2)
    elif cornerUL < R:
        theta = abs( abs(math.asin( np.clip(abs(cirY - y1) / R,-1.0,1.0) )) - abs(math.asin ( np.clip(abs(cirY - y2) / R,-1.0,1.0) )) )
        t2 = math.pi/2 -  abs( abs(math.asin( np.clip(abs(cirY - y2) / R,-1.0,1.0) )) + abs(math.asin ( np.clip(abs(cirX - x2) / R,-1.0,1.0) )) )
#DEBUG  #print "Theta: %f 90-t2: %f" % (theta,t2)
#DEBUG  #print "Corner 3: %f Theta: %f x2: %f y2: %f"% (cornerUL,theta,x2,y2)
    elif cornerLL < R:
        theta = abs( abs(math.asin( np.clip(abs(cirX - x1) / R,-1.0,1.0) )) - abs(math.acos ( np.clip(abs(cirY - y1) / R,-1.0,1.0) )) )
        t2 = math.pi/2 -  abs( abs(math.asin( np.clip(abs(cirY - y1) / R,-1.0,1.0) )) + abs(math.asin ( np.clip(abs(cirX - x1) / R,-1.0,1.0) )) )
        #print "Theta: %f 90-t2: %f" % (theta,t2)

    A2 = abs( (x2-x1)*(y2-y1) / 2 )
    A1 = R**2 * abs( theta/2 - (math.cos(theta/2) * math.sin(theta/2)) )
    A = A1 + A2
    if A > 1:
        raise "Onecorner A > 1"
    elif A < 0:
        raise "Onecorner A < 0"
    return A

def threecorner(innerX,innerY,outerX,outerY,R,cirX,cirY):
    if innerY < outerY: 
        y3 = cirY + np.sqrt( R**2 - (cirX-outerX)**2 )
    elif innerY > outerY: 
        y3 = cirY - np.sqrt( R**2 - (cirX-outerX)**2 )
    if innerX < outerX:
        x3 = cirX + np.sqrt( R**2 - (cirY-outerY)**2 )
    elif innerX > outerX:
        x3 = cirX - np.sqrt( R**2 - (cirY-outerY)**2 )
    A2 = abs( (x3-outerX)*(y3-outerY) / 2 )
    A3 = abs( (outerX-innerX)*(outerY-innerY) ) - 2*A2 
    theta = abs( math.acos( abs(cirX - outerX) / R ) - math.acos ( abs(cirX - x3) / R ) )
    A1 = R**2 * abs( theta/2 - (math.cos(theta/2) * math.sin(theta/2)) )
    A = A1 + A2 + A3
#DEBUGp rint "Three corner: A1:%f A2:%f A3:%f theta:%f" % (A1,A2,A3,theta)
    if A > 1:
        raise "Threecorner A > 1"
    elif A < 0:
        raise "Threecorner A < 0"
    return A

def twocorner(x1,x2,y1,y2,cornerUR,cornerLR,cornerUL,cornerLL,R,cirX,cirY):
    if (cornerUR < R and cornerLR < R):
        x4 = cirX - np.sqrt( R**2 - (cirY - y1)** 2 )
        x3 = cirX - np.sqrt( R**2 - (cirY - y2)** 2 )
        A2 = abs( (x3-x4)*(y1-y2) / 2 )
        if abs(x2-x3)<abs(x2-x4):
            A3 = abs( (x2-x3)*(y1-y2) )
        else:
            A3 = abs( (x2-x4)*(y1-y2) )
        theta = abs( math.asin( np.clip((cirY - y1) / R,-1.0,1.0) ) - math.asin ( np.clip((cirY - y2) / R,-1.0,1.0) ) )
        A1 = R**2 * abs( theta/2 - abs(math.cos(theta/2) * math.sin(theta/2)) )
    elif cornerLR < R and cornerLL < R:
        y4 = cirY + np.sqrt( R**2 - (cirX - x1)** 2 )
        y3 = cirY + np.sqrt( R**2 - (cirX - x2)** 2 )
        A2 = abs( (y3-y4)*(x1-x2) / 2 )
        if abs(y1-y3)<abs(y1-y4):
            A3 = abs( (y1-y3)*(x1-x2) )
        else:
            A3 = abs( (y1-y4)*(x1-x2) )
        theta = abs( abs( math.asin( np.clip((cirX - x1) / R,-1.0,1.0) ) ) - abs ( math.asin ( np.clip((cirX - x2) / R,-1.0,1.0) ) ) )
        A1 = R**2 * abs( theta/2 - abs(math.cos(theta/2) * math.sin(theta/2)) )
    elif cornerUL < R and cornerLL < R:
        x4 = cirX + np.sqrt( R**2 - (cirY - y1)** 2 )
        x3 = cirX + np.sqrt( R**2 - (cirY - y2)** 2 )
        A2 = abs( (x3-x4)*(y1-y2) / 2 )
        if abs(x1-x3)<abs(x1-x4):
            A3 = abs( (x1-x3)*(y1-y2) )
        else:
            A3 = abs( (x1-x4)*(y1-y2) )
        theta = abs( math.asin( np.clip((cirY - y1) / R,-1.0,1.0) ) - math.asin ( np.clip((cirY - y2) / R,-1.0,1.0) ) )
        A1 = R**2 * abs( theta/2 - abs(math.cos(theta/2) * math.sin(theta/2)) )
    elif cornerUL < R and cornerUR < R:
        y4 = cirY - np.sqrt( R**2 - (cirX - x1)** 2 )
        y3 = cirY - np.sqrt( R**2 - (cirX - x2)** 2 )
        A2 = abs( (y3-y4)*(x1-x2) / 2 )
        if abs(y2-y3)<abs(y2-y4):
            A3 = abs( (y2-y3)*(x1-x2) )
        else:
            A3 = abs( (y2-y4)*(x1-x2) )
        theta = abs( math.asin( np.clip((cirX - x1) / R,-1.0,1.0) ) - math.asin ( np.clip((cirX - x2) / R,-1.0,1.0) ) )
        A1 = R**2 * abs( theta/2 - abs(math.cos(theta/2) * math.sin(theta/2)) )
    else:
        print "Corners 1 %f 2 %f 3 %f 4 %f" % (cornerUR,cornerLR,cornerUL,cornerLL)
        raise ValueError("Error x1 %2.2f x2 %2.2f y1 %2.2f y2 %2.2f radius %2.3f" % (x1,x2,y1,y2,R))
    A = A1 + A2 + A3
    if A > 1:
        #WARNING: THIS IS A HACK!
        return 1
        print "Input parameters: " 
        print "x1,x2,y1,y2,cornerUR,cornerLR,cornerUL,cornerLL,R,cirX,cirY"
        print x1,x2,y1,y2,cornerUR,cornerLR,cornerUL,cornerLL,R,cirX,cirY
        print "A1 %f A2 %f A3 %f  theta: %f" % (A1,A2,A3,theta)
        raise ValueError("Twocorner A > 1")
    elif A < 0:
        print "A1 %f A2 %f A3 %f" % (A1,A2,A3)
        raise "Twocorner A < 0"
    elif theta <= 0:
        raise "Theta <= 0"
    return A

def circmask(cx,cy,r,xsize,ysize):
#DEBUG  print cx,cy,r,xsize,ysize
    def f(x1,y1):
        return float(pixelarea(float(x1),float(y1),float(cx),float(cy),float(r)))
    f2 = np.vectorize(f)
    return np.fromfunction(f2,(xsize,ysize))

