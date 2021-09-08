# Functions for generating VRML scenes
# The *.wrl functions are internal, not meant to be called from outside.

truename <- function(f) {
  # turn a path into an absolute path
  if((substring(f,2,2) == ":") || (substring(f,1,1) %in% c("\\","/")))
    return(f)
  paste(getwd(),f,sep="\\")
}

vrml.new.name <- function() {
  assign("vrml.counter",vrml.counter + 1,globalenv())
  paste("obj",vrml.counter,sep="")
}

vrml.color <- function(col) {
  t(col2rgb(col))/255
}
vrml.material <- list(diffuseColor=c(1,1,1),
                      emissiveColor=c(0,0,0),
                      specularColor=c(0,0,0),
                      transparency=0,shininess=0.2,
                      ambientIntensity=0.2)
vrml.material.wrl <- function(diffuseColor=vrml.material$diffuseColor,
                              emissiveColor=vrml.material$emissiveColor,
                              transparency=vrml.material$transparency,
                              specularColor=vrml.material$specularColor,
                              shininess=vrml.material$shininess,
                              ambientIntensity=vrml.material$ambientIntensity,
                              ...,file=vrml.file) {
  # writes an Appearance node
  if(is.null(diffuseColor)) diffuseColor = vrml.material$diffuseColor
  cat("Appearance { material Material {",file=file)
  if(any(emissiveColor != 0)) cat(" emissiveColor",emissiveColor,file=file)
  if(any(specularColor != 0)) cat(" specularColor",specularColor,file=file)
  if(any(diffuseColor != 0.8)) cat(" diffuseColor",diffuseColor,file=file)
  if(transparency != 0) cat(" transparency",transparency,file=file)
  if(ambientIntensity != 0.2) cat(" ambientIntensity",ambientIntensity,file=file)
  # shininess does nothing without specularColor
  if(shininess != 0.2) cat(" shininess",shininess,file=file)
  cat(" } } ",file=file)
}
vrml.clear.material.cache <- function() {
  assign("vrml.materials",list(),globalenv())
}
vrml.as.material <- function(diffuseColor=vrml.material$diffuseColor,
                             emissiveColor=vrml.material$emissiveColor,
                             transparency=vrml.material$transparency,
                             specularColor=vrml.material$specularColor,
                             shininess=vrml.material$shininess,
                             ambientIntensity=vrml.material$ambientIntensity,
                             ...) {
  if(is.null(diffuseColor)) diffuseColor = vrml.material$diffuseColor
  list(diffuseColor=diffuseColor,
       emissiveColor=emissiveColor,
       transparency=transparency,
       specularColor=specularColor,
       shininess=shininess,
       ambientIntensity=ambientIntensity)
}
vrml.material.key <- function(...) {
  color.spec = vrml.as.material(...)
  color.spec = lapply(color.spec, function(x) {attributes(x) = NULL; x})
  paste(deparse(color.spec),collapse=" ")
}
vrml.material.wrl.cache <- function(...,file=vrml.file) {
  # writes an appearance field, with USE or DEF as appropriate
  if(!exists("vrml.materials")) vrml.clear.material.cache()
  key = vrml.material.key(...)
  the.appearance = vrml.materials[[key]]
  if(is.null(the.appearance)) {
    the.appearance = vrml.new.name()
    cat("appearance DEF",the.appearance,"",file=file)
    vrml.material.wrl(...,file=file)
    vrml.materials[[key]] = the.appearance
    assign("vrml.materials",vrml.materials,globalenv())
  } else {
    cat("appearance USE",the.appearance,"",file=file)
  }
}
vrml.material.keys <- function(color,...) {
  key = c()
  for(i in 1:nrow(color)) {
    key[i] = vrml.material.key(color[i,],...)
  }
  key
}

vrml.set.scale <- function(xlim,ylim,zlim,scale=c(1,1,1)) {
  # sets the globals vrml.offset and vrml.scale,
  # each vectors of 3 numbers
  s = scale/c(diff(range(xlim)),diff(range(ylim)),diff(range(zlim)))
  offset = -c(xlim[1],ylim[1],zlim[1])
  assign("vrml.scale",s,globalenv())
  assign("vrml.offset",offset,globalenv())
  assign("vrml.symbol.radius",0.01,globalenv())
}
vrml.to.scale <- function(p) {
  # p is matrix of three columns
  if(ncol(p) != 3) stop("p must be a matrix of 3 columns")
  p = p + rep.row(vrml.offset,nrow(p))
  scale.cols(p, vrml.scale)
}

# http://www.cs.cmu.edu/afs/cs/academic/class/16741-s02/www/lecture6.pdf
# http://www.cs.cmu.edu/afs/cs/academic/class/16741-s02/www/lecture7.pdf
vrml.quaternion.rotation <- function(rot) {
  rot[1:3] = rot[1:3]/sqrt(sum(rot[1:3]^2))
  c(cos(rot[4]/2),sin(rot[4]/2)*rot[1:3])
}
vrml.rotation.quaternion <- function(q) {
  qv = q[2:4]
  if(all(qv == 0)) c(0,1,0,0)
  else c(qv/norm(qv), 2*atan2(norm(qv),q[1]))
}
vrml.compose.quaternion <- function(q1,q2) {
  q1v = q1[2:4]
  q2v = q2[2:4]
  c(q1[1]*q2[1] - sum(q1v*q2v),q1[1]*q2v+q2[1]*q1v+vec.cross(q1v,q2v))
}
vrml.compose.rotation <- function(rot1,rot2) {
  q1 = vrml.quaternion.rotation(rot1)
  q2 = vrml.quaternion.rotation(rot2)
  q = vrml.compose.quaternion(q1,q2)
  vrml.rotation.quaternion(q)
}
vec.cross <- function(p1,p2) {
  # cross-product of vectors
  # http://mathworld.wolfram.com/CrossProduct.html
  cbind(p1[2]*p2[3] - p1[3]*p2[2],
        p1[3]*p2[1] - p1[1]*p2[3],
        p1[1]*p2[2] - p1[2]*p2[1])
}
vec.angle <- function(p1,p2) {
  # angle between two vectors
  acos(sum(p1*p2)/sqrt(sum(p1*p1)*sum(p2*p2)))
}
vrml.rotate <- function(p,rot) {
  angle = rot[4]
  v = rot[1:3]
  vp = vec.cross(v,p)
  vvp = vec.cross(v,vp)
  p + sin(angle)*vp + (1-cos(angle))*vvp
}
vrml.rotation.between <- function(p1,p2) {
  # returns a four-element vector specifying an axis-angle rotation
  # that takes p1 to p2
  if(!is.null(dim(p1)) && nrow(p1) > 1) {
    rot = vrml.rotation.between(p1[1,,drop=F],p2[1,,drop=F])
    p1[2,] = vrml.rotate(p1[2,],rot)
    rot2 = vrml.rotation.between(p1[2,,drop=F],p2[2,,drop=F])
    vrml.compose.rotation(rot2,rot)
  } else {
    if(all(p1 == 0) || all(p2 == 0)) stop("impossible rotation (to/from the origin)")
    axis = vec.cross(p1,p2)
    if(sum(axis*axis) == 0) axis = cbind(0,1,0)
    else axis = axis/sqrt(sum(axis*axis))
    c(axis,vec.angle(p1,p2))
  }
}

vrml.viewpoint.wrl <- function(pos,lookat,file=vrml.file) {
  default.orient = cbind(0,0,-1)
  rot = vrml.rotation.between(default.orient,lookat-pos)
  cat("Viewpoint { position",pos,"orientation",rot,"}\n",file=file)
}

vrml.file = ""

vrml.open <- function(xlim,ylim,zlim,scale=c(1,1,1),light=F,bg=c(1,1,1),...,
                      file.name=NULL) {
  if(is.null(file.name)) file.name = file.path(tempdir(),"Rplot")
  if(substring(file.name,nchar(file.name)-3,nchar(file.name)) != ".wrl") {
    file.name = paste(file.name,"wrl",sep=".")
  }
  if(nchar(file.name) > 4) {
    con = file(file.name,"w")
  } else {
    con = ""
  }
  assign("vrml.file",con,globalenv())
  assign("vrml.file.name",file.name,globalenv())
  assign("vrml.counter",0,globalenv())
  vrml.clear.material.cache()
  cat("#VRML V2.0 utf8\n",file=vrml.file)
  #cat("WorldInfo { title \"",match.call()[[1]],"\" }\n",sep="",file=vrml.file)
  cat("NavigationInfo { type \"EXAMINE\" speed 10 ",file=vrml.file)
  if(light) cat("headlight FALSE ",file=vrml.file)
  cat("}\n",file=vrml.file)
  # white background is easier on the eyes
  if(length(bg) != 3) bg = col2rgb(bg)[,]
  cat("Background { skyColor [",bg,"] }\n",file=vrml.file)
  cat("DEF black Appearance { material Material {} }\n",file=vrml.file)
  # camera looks at the center of the back wall
  p = cbind(0.5,-1,0.5)*scale - cbind(0,0.5,0)
  lookat = cbind(0.5,0,0.5)*scale
  vrml.viewpoint.wrl(p,lookat)
  vrml.set.scale(xlim,ylim,zlim,scale)
}

vrml.close <- function() {
  if(inherits(vrml.file,"connection")) {
    close(vrml.file)
    cat("wrote to",vrml.file.name,"\n")
    shell.exec = if(exists("shell.exec")) shell.exec else browseURL
    shell.exec(truename(vrml.file.name))
  }
}

# baseline specifies the plane in which text is drawn
vrml.baseline = rbind(c(1,0,0),c(0,1,0))
vrml.text.wrl <- function(p,labels,adj=c(0,0),color=NULL,
                          baseline=vrml.baseline,cex=1,billboard=F,...,
                          file=vrml.file) {
  # if only one vector given for baseline, use the default height vector
  if(nrow(baseline) == 1) baseline = rbind(baseline,vrml.baseline[2,])
  rot = vrml.rotation.between(vrml.baseline,baseline)
  s = cex*vrml.symbol.radius*6
  if(adj[1] > 0) {
    # BUG: this requires an existing plot
    w = strwidth(labels,units="i")*5
    dim(w) = c(length(w),1)
    p = p - (w %*% baseline[1,,drop=F]*adj[1]*s)
  }
  if(adj[2] > 0) {
    h = strheight(labels,units="i")*6
    dim(h) = c(length(h),1)
    p = p - (h %*% baseline[2,,drop=F]*adj[2]*s)
  }
  for(i in 1:nrow(p)) {
    cat("Transform { translation",p[i,],"\n",
        "rotation",rot,"\n",
        "scale",rep(s,3),"\n",
        "children ",file=file)
    if(billboard) {
      cat("Billboard { axisOfRotation 0 0 0 children ",file=file)
    }
    cat("Shape { ",file=file)
    color.i = if(is.null(color) || nrow(color) == 1) color else color[i,]
    vrml.material.wrl.cache(color.i,...,file=file)
    txt = paste("\"",labels[i],"\"",sep="")
    cat("geometry Text { string",txt,"} }",file=file)
    if(billboard) cat(" } ",file=file)
    cat("}\n",file=file)
  }
}
vrml.text <- function(x,y,z,col=1,...) {
  p = vrml.to.scale(cbind(x,y,z))
  vrml.text.wrl(p,color=vrml.color(col),...)
}

vrml.sphere <- function(x,y,z,cex=1,col=1,...) {
  p = vrml.to.scale(cbind(x,y,z))
  r = cex*vrml.symbol.radius
  vrml.sphere.wrl(p,r,color=vrml.color(col),...)
}
vrml.sphere.wrl <- function(p,r,color=NULL,...,
                            file=vrml.file) {
  if(length(r) == 1) r = rep(r,nrow(p))
  for(i in 1:nrow(p)) {
    cat("Transform { translation",p[i,],"\n",
        "scale",rep(r[i],3),"\n",
        "children Shape { ",file=file)
    color.i = if(is.null(color) || nrow(color) == 1) color else color[i,]
    vrml.material.wrl.cache(color.i,...,file=file)
    if(i == 1) {
      the.shape = vrml.new.name()
      cat("geometry DEF",the.shape,"Sphere {} ",file=file)
    } else {
      cat("geometry USE",the.shape,"",file=file)
    }
    cat("}\n}\n",file=file)
  }
}
vrml.cube <- function(x,y,z,s=c(1,1,1),cex=1,col=1,...) {
  p = vrml.to.scale(cbind(x,y,z))
  s = s*vrml.symbol.radius*cex
  vrml.cube.wrl(p,s,color=vrml.color(col),...)
}
vrml.cube.wrl <- function(p,s=c(1,1,1),color=NULL,...,
                          file=vrml.file) {
  # this would be more efficient with an IndexedFaceSet
  #if(is.null(dim(color))) dim(color) = c(1,3)
  if(length(s) == 1) s = rep(s,3)
  p = na.omit(p)
  omitted = attr(p,"na.action")
  if(is.null(dim(s)) || nrow(s)==1) s = rep.row(s,nrow(p))
  else if(!is.null(omitted)) s = s[omitted,]
  for(i in 1:nrow(p)) {
    cat("Transform { translation",p[i,],"\n",
        "scale",s[i,],"\n",
        "children Shape { ",file=file)
    color.i = if(is.null(color) || nrow(color) == 1) color else color[i,]
    vrml.material.wrl.cache(color.i,...,file=file)
    if(i == 1) {
      the.shape = vrml.new.name()
      cat("geometry DEF",the.shape,"Box {} ",file=file)
    } else {
      cat("geometry USE",the.shape,"",file=file)
    }
    cat("}\n}\n",file=file)
  }
}

vrml.cylinder.wrl <- function(p1,p2,r=vrml.symbol.radius,file=vrml.file) {
  midp = (p1+p2)/2
  dp = p2-p1
  h = sqrt(sum(dp*dp))
  rot = vrml.rotation.between(cbind(0,1,0),dp)
  cat("Transform { translation",midp[1,],"\n",
      "rotation",rot,"\n",
      "children Shape { geometry Cylinder { radius",r,"height",h,"} }\n",
      "}\n",
      file=file)
}

vrml.axes <- function(xlim,ylim,zlim,xlab=NULL,ylab=NULL,zlab=NULL,cex=2,...,
                      file=vrml.file) {
  # origin
  po = vrml.to.scale(cbind(xlim[1],ylim[1],zlim[1]))
  vrml.cube.wrl(po,vrml.symbol.radius,...,file=file)
  px = vrml.to.scale(cbind(xlim[2],ylim[1],zlim[1]))
  py = vrml.to.scale(cbind(xlim[1],ylim[2],zlim[1]))
  pz = vrml.to.scale(cbind(xlim[1],ylim[1],zlim[2]))
  if(F) {
    # cylinders for axes
    vrml.cylinder.wrl(po,px,file=file)
    vrml.cylinder.wrl(po,py,file=file)
    vrml.cylinder.wrl(po,pz,file=file)
  } else {
    vrml.line.wrl(po,px,file=file)
    vrml.line.wrl(po,py,file=file)
    vrml.line.wrl(po,pz,file=file)
  }
  # axis names
  mar = vrml.symbol.radius
  if(!is.null(xlab)) {
    vrml.text.wrl(px/2-pz*mar,xlab,cex=cex,adj=c(0.5,1),base=rbind(px,pz),file=file)
    #vrml.text.wrl(px/2-py*mar+pz,xlab,cex=cex,adj=c(0.5,1),base=px,file=file)
    #vrml.text.wrl(px/2-py*mar,xlab,cex=cex,adj=c(0.5,1),base=-px,file=file)
  }
  if(!is.null(ylab)) {
    vrml.text.wrl(py/2-pz*mar,ylab,cex=cex,adj=c(0.5,1),base=rbind(py,pz),file=file)
    #vrml.text.wrl(py/2-px*mar+pz,ylab,cex=cex,adj=c(0.5,1),base=py,file=file)
    #vrml.text.wrl(py/2-px*mar,ylab,cex=cex,adj=c(0.5,1),base=-py,file=file)
  }
  if(!is.null(zlab)) {
    vrml.text.wrl(pz/2-px*mar,zlab,cex=cex,adj=c(0.5,0),base=rbind(pz,-px),file=file)
    #vrml.text.wrl(pz/2-py*mar+px,zlab,cex=cex,adj=c(0.5,1),base=-pz,file=file)
    #vrml.text.wrl(pz/2-py*mar,zlab,cex=cex,adj=c(0.5,1),base=pz,file=file)
  }
}
vrml.line.wrl <- function(p1,p2,file=vrml.file) {
  cat("Shape { appearance USE black geometry IndexedLineSet {\n",
      "coord Coordinate {\n",
      "point [",t(rbind(p1,p2)),"]\n }\n",
      "coordIndex [",c(0,1,-1),"]\n}}\n",
      file=file)
}
vrml.points <- function(x,y,z,col=1,...) {
  p = vrml.to.scale(cbind(x,y,z))
  vrml.points.wrl(p,color=vrml.color(col),...)
}
vrml.points.wrl <- function(p,color=NULL,...,
                            file=vrml.file) {
  cat("Shape { ",file=file)
  color.i = if(is.null(color) || nrow(color) == 1) color else color[1,]
  vrml.material.wrl.cache(color.i,...,file=file)
  cat("geometry PointSet {\n",file=file)
  if(nrow(color) == nrow(p)) {
    cat("color Color { color [",t(color),"] }\n",file=file)
  }
  cat("coord Coordinate { point [",t(p),"] }\n",
      "}}\n",
      file=file)
}

vrml.box.wrl <- function(p,file=vrml.file) {
  points = array(0,c(8,3))
  g = data.matrix(expand.grid(1:2,1:2,1:2))
  for(i in 1:nrow(g)) {
    for(j in 1:3) points[i,j] = p[g[i,j],j]
  }
  index = c(0,1,5,5,4,0,-1,2,3,7,6,2,-1,
            0,2,-1,1,3,-1,5,7,-1,4,6,-1)
  cat("Shape { appearance USE black geometry IndexedLineSet {\n",
      "coord Coordinate {\n",
      "point [",t(points),"]\n }\n",
      "coordIndex [",index,"]\n}}\n",
      file=file)
}
vrml.box <- function(xlim,ylim,zlim,...) {
  vrml.box.wrl(vrml.to.scale(cbind(xlim,ylim,zlim)),...)
}

#' View points in 3D
#'
#' Creates a VRML scene with the given points.
#' @aliases vrml.plot3d vrml.plot3d.default vrml.plot3d.data.frame vrml.plot3d.formula
#' @param x,y,z numeric vectors giving the coordinates of the points.
#' @param frame a data frame with three columns.
#' @param data an environment in which to lookup x,y,z.
#' @param xlab,ylab,zlab}{axis labels.
#' @param pch a numeric code for the graphic object to place at
#'     each location.  See details.
#' @param col a single color or vector of colors, used to color the points.
#' @param scale a vector of 3 numbers, defining the size of the box into
#'     which the points are placed.  You can control the aspect ratio this way.
#' @param file.name a filename (with or without the .wrl extension) to
#'     receive the VRML.  If \code{NULL}, a temporary name is chosen.
#' @param cex a scale factor for the objects.
#' @param cex.axis a scale factor for the axis labels.
#' @param light If \code{TRUE}, the scene will contain a light at each
#'     corner of the cube, with no headlight.
#'     Otherwise, the scene will contain no lights and rely on the
#'     headlight only.
#' @param axes If \code{TRUE}, a box with axis labels is drawn around the
#'     points.
#' @param ... additional arguments to pass to the internal drawing routines,
#'     or to pass through to \code{vrml.plot3d.default}
#'     (for the \code{data.frame} and \code{formula} methods).
#' @details
#'   A VRML scene description file is created, and opened with a browser if
#'   one is available.
#'
#'   VRML is a standard language for describing 3D scenes, and the file
#'   produced by this function should be portable across all machines
#'   with a VRML browser.
#'
#'   If \code{pch=0}, a cube is placed at each location.
#'   If \code{pch=1}, a sphere is placed at each location.
#'   If \code{pch="."}, a dot is placed at each location.
#'   The latter option is the fastest, and will be chosen automatically if
#'   the scale is small or the dataset is large.
#' @author Tom Minka
#' @references
#'   \url{http://web3d.vapourtech.com/tutorials/vrml97/}
#'
#'   \url{www.wed3d.org/resources/vrml_ref_manual/Book.html}
#' @examples
#' data(Housing)
#' w = pca(HousingT,k=3)
#' x = project(HousingT,w)
#' plot.new()
#' vrml.plot3d(x)
#' @export
vrml.plot3d <- function(object, ...) UseMethod("vrml.plot3d")


#' @export
vrml.plot3d.formula <- function(formula,data=parent.frame(),...) {
  x = model.frame.default(formula,data,na.action=na.pass)
  vrml.plot3d.data.frame(x,...)
}


#' @export
vrml.plot3d.data.frame <- function(x,labels=NULL,...,xlab,ylab,zlab) {
  resp = response.var(x)
  pred = predictor.vars(x)
  if(missing(xlab)) xlab = pred[1]
  if(missing(ylab)) ylab = pred[2]
  if(missing(zlab)) zlab = resp
  if(identical(labels,TRUE)) labels = rownames(x)
  vrml.plot3d.default(x[,pred[1]],x[,pred[2]],x[,resp],labels=labels,
                      xlab=xlab,ylab=ylab,zlab=zlab,...)
}


#' @export
vrml.plot3d.default <- function(x,y,z,data=parent.frame(),labels=NULL,
                                xlab,ylab,zlab,
                                pch=0,col=3,
                                scale=c(1,1,1),asp=NULL,
                                file.name=NULL,
                                cex=1,cex.axis=2,light=T,axes=T,type="p",...) {
  if(missing(xlab)) xlab <- deparse(substitute(x))
  if(missing(ylab)) ylab <- deparse(substitute(y))
  if(missing(zlab)) zlab <- deparse(substitute(z))
  x <- eval(substitute(x),data)
  y <- eval(substitute(y),data)
  z <- eval(substitute(z),data)

  xlim = range(x,na.rm=T); ylim = range(y,na.rm=T); zlim = range(z,na.rm=T)
  xlim = extend.pct(xlim); ylim = extend.pct(ylim); zlim = extend.pct(zlim)
  if(!is.null(asp)) {
    if(length(asp) == 1) asp = rep(asp,2)
    scale = c(diff(xlim),diff(ylim)*asp[1],diff(zlim)*asp[2])
    scale = scale/diff(xlim)
  }
  vrml.open(xlim,ylim,zlim,scale,light,file.name=file.name,...)
  if(light) {
    cat("PointLight { ambientIntensity 0.5 intensity 0.25 location 0 0 1 }\n",
        file=vrml.file)
    cat("PointLight { intensity 0.25 location 1 1 1 }\n",
        file=vrml.file)
    cat("PointLight { intensity 0.25 location 0 1 1 }\n",
        file=vrml.file)
    cat("PointLight { intensity 0.25 location 1 0 1 }\n",
        file=vrml.file)
  }
  if(axes) {
    cat("Group { children [\n",file=vrml.file)
    vrml.box(xlim,ylim,zlim)
    vrml.axes(xlim,ylim,zlim,xlab,ylab,zlab,cex=cex.axis,...)
    cat("]}\n",file=vrml.file)
  }
  cat("Group { children [\n",file=vrml.file)
  if(any(scale/cex >= 5) || length(x) > 5000) {
    pch = "."
    cat("Using pch=. to save time\n")
  }
  if(is.null(labels)) {
    if(pch == ".") {
      vrml.points(x,y,z,col=col,...)
    } else if(pch==1) {
      vrml.sphere(x,y,z,col=col,cex=cex,...)
    } else {
      # pch=0 is a cube
      vrml.cube(x,y,z,col=col,cex=cex,...)
    }
  } else {
    vrml.text(x,y,z,labels=labels,col=col,cex=cex,billboard=T,...)
  }
  if(type == "o") {
    p = vrml.to.scale(cbind(x,y,z))
    cat("Shape { appearance USE black geometry IndexedLineSet {\n",
        "coord Coordinate {\n",
        "point [",t(p),"]\n }\n",
        "coordIndex [",c(1:nrow(p),0)-1,"]\n}}\n",
        file=vrml.file)
  }
  cat("]}\n",file=vrml.file)
  vrml.close()
}


vrml.color.plot <- function(x,y,nlevels=4,color.palette=YlGnBu.colors,...) {
  if(missing(y)) {
    resp = response.var(x)
    y = x[[resp]]
    pred = predictor.terms(x)
    x = x[pred]
  }
  if(!is.factor(y)) y = cut.quantile(y,nlevels)
  else {
    color.palette = default.colors
    nlevels = length(levels(y))
  }
  if(is.function(color.palette)) color.palette = color.palette(nlevels)
  col = color.palette[y]
  vrml.plot3d(formula(x),x,col=col,light=F,...)
}

#' View a surface in 3D
#'
#' Creates a VRML scene with a shaded surface.
#' @aliases vrml.surface vrml.surface.default vrml.surface.loess
#' @param x,y locations of grid lines at which the values in \code{z} are
#'     measured.
#' @param z}{a matrix containing the values to be plotted (\code{NA}s are
#'                                                          allowed).
#' @param xlab,ylab,zlab}{axis labels.  If \code{NULL}, taken from the
#'     deparsed expressions for \code{x,y,z}.
#' @param col the color of the surface.
#' @param scale a vector of 3 numbers, defining the size of the box into
#'     which the surface is placed.  You can control the aspect ratio this way.
#' @param file.name a filename (with or without the .wrl extension) to
#'     receive the VRML.  If \code{NULL}, a temporary name is chosen.
#' @param cex.axis a scale factor for the axis labels.
#' @param light If \code{TRUE}, the scene will contain its own light
#'     and no headlight.
#'     Otherwise, the scene will contain no lights and rely on the
#'     headlight only.
#' @param border If \code{TRUE}, a grid will be drawn on top of the surface.
#' @param creaseAngle a parameter controlling the smoothness of the
#'     surface.  When it is small, faces at large angles to each other will
#'     create visible "creases".
#' @param ... additional arguments for internal drawing routines.
#' @details
#'   A VRML scene description file is created, and opened with a browser if
#'   one is available.
#'
#'   VRML is a standard language for describing 3D scenes, and the file
#'   produced by this function should be portable across all machines
#'   with a VRML browser.
#'
#'   This function is similar to \code{\link{persp}} except it creates a 3D
#'   scene which can be manipulated with a VRML viewer.
#' @author Tom Minka
#' @seealso \code{\link{persp}},\code{\link{vrml.plot3d}}
#' @examples
#' #  The Obligatory Mathematical surface.
#' #  Rotated sinc function.
#' x <- seq(-10, 10, length= 30)
#' y <- x
#' f <- function(x,y) { r <- sqrt(x^2+y^2); 10 * sin(r)/r }
#' z <- outer(x, y, f)
#' z[is.na(z)] <- 1
#' vrml.surface(x,y,z)
#' vrml.surface(x,y,z,border=T)
#' @export
vrml.surface <- function(object,...) UseMethod("vrml.surface")

#' @export
vrml.surface.loess <- function(object,res=20,
                               xlim,ylim,zlim,clip=T,xlab,ylab,zlab,...) {
  if(length(res) == 1) res <- rep(res,2)
  x <- model.frame(object)
  resp <- response.var(object)
  pred <- predictor.vars(object)
  if(missing(xlim)) xlim <- range(x[[pred[1]]],na.rm=T)
  x1 <- seq(xlim[1],xlim[2],length=res[1])
  if(missing(ylim)) ylim <- range(x[[pred[2]]],na.rm=T)
  x2 <- seq(ylim[1],ylim[2],length=res[2])
  xt <- expand.grid(x1,x2)
  names(xt) <- pred[1:2]
  z <- predict(object,xt)
  if(!missing(zlim)) {
    z[z < zlim[1]] <- zlim[1]
    z[z > zlim[2]] <- zlim[2]
  }
  if(clip != F) {
    if(clip == T) {
      i <- chull(x[pred])
      clip <- x[pred][i,]
    }
    z[!in.polygon(clip,xt)] <- NA
  }
  dim(z) <- c(length(x1),length(x2))
  if(missing(xlab)) xlab <- pred[1]
  if(missing(ylab)) ylab <- pred[2]
  if(missing(zlab)) zlab <- resp
  vrml.surface.default(x1,x2,z,xlab=xlab,ylab=ylab,zlab=zlab,...)
}

#' @export
vrml.surface.default <- function(x,y,z,xlab=NULL,ylab=NULL,zlab=NULL,
                                 col="gray",scale=c(1,1,1),file.name=NULL,
                                 cex.axis=2,light=F,...) {
  if(is.null(xlab)) xlab <- deparse(substitute(x))
  if(is.null(ylab)) ylab <- deparse(substitute(y))
  if(is.null(zlab)) zlab <- deparse(substitute(z))
  #x <- eval(substitute(x),data)
  #y <- eval(substitute(y),data)
  #z <- eval(substitute(z),data)

  xlim = range(x,na.rm=T); ylim = range(y,na.rm=T); zlim = range(z,na.rm=T)
  xlim = extend.pct(xlim); ylim = extend.pct(ylim); zlim = extend.pct(zlim)
  if(identical(scale,F)) {
    scale = c(diff(range(xlim)),diff(range(ylim)),diff(range(zlim)))
    scale = scale/sum(scale)*3
  }
  vrml.open(xlim,ylim,zlim,scale,light,file.name=file.name,...)
  if(light) {
    # dir 1 1 -1
    cat("DirectionalLight { ambientIntensity 0 intensity 1",
        "direction 1 1 -0.1 }\n",file=vrml.file)
  }
  cat("Group { children [\n",file=vrml.file)
  vrml.box(xlim,ylim,zlim)
  vrml.axes(xlim,ylim,zlim,xlab,ylab,zlab,cex=cex.axis,...)
  cat("]}\n",file=vrml.file)
  cat("Group { children [\n",file=vrml.file)
  vrml.surface.wrl(x,y,z,color=vrml.color(col),...)
  cat("]}\n",file=vrml.file)
  vrml.close()
}

#' @export
vrml.surface.wrl <- function(x,y,z,border=F,
                             creaseAngle=10,...,file=vrml.file) {
  # turn off LOD in the viewer or you will get weird effects
  # could also use an ElevationGrid
  points = expand.grid(x,y)
  points = cbind(points,as.vector(z))
  points = vrml.to.scale(points)
  face = 0
  nx = length(x)
  ny = length(y)
  index = array(0,c((nx-1)*(ny-1)*2,4))
  for(i in 0:(ny-2)) {
    for(j in 0:(nx-2)) {
      face = face + 1
      index[face,] = c(i*nx + j+1, i*nx + j, (i+1)*nx + j, -1)
      if(any(is.na(points[index[face,1:3]+1,3]))) face = face - 1
      face = face + 1
      index[face,] = c(i*nx + j+1, (i+1)*nx + j, (i+1)*nx + j+1, -1)
      if(any(is.na(points[index[face,1:3]+1,3]))) face = face - 1
    }
  }
  index = index[1:face,]
  cat("Shape { appearance ",file=file)
  vrml.material.wrl(...,file=file)
  cat("geometry IndexedFaceSet { solid FALSE creaseAngle",creaseAngle,"\n",
      "coord Coordinate {\n",
      "point [",t(na.dummy(points)),"]\n }\n",
      "coordIndex [",t(index),"]\n}}\n",file=file)

  if(border) {
    points[,3] = points[,3] + 4e-3
    index = array(0,dim=c((nx-1)*(ny-1),5))
    face = 0
    for(i in 0:(ny-2)) {
      for(j in 0:(nx-2)) {
        face = face + 1
        index[face,] = c(i*nx + j+1, i*nx + j, (i+1)*nx + j, (i+1)*nx + j+1, -1)
        if(any(is.na(points[index[face,1:4]+1,3]))) face = face - 1
      }
    }
    index = index[1:face,]
    cat("Shape { appearance ",file=file)
    # only emissiveColor is used
    vrml.material.wrl(emissiveColor=vrml.color(1),file=file)
    cat("geometry IndexedLineSet {",
        "coord Coordinate {\n",
        "point [",t(na.dummy(points)),"]\n }\n",
        "coordIndex [",t(index),"]\n}}\n",file=file)
  }
}


vrml.array.wrl <- function(m,file=vrml.file) {
  if(is.data.frame(m)) m = data.matrix(m)
  for(i in 1:nrow(m)) {
    cat(m[i,],", ",file=file)
  }
}
vrml.cone.grid <- function() {
  ntheta = 24
  points = array(0,dim=c(6*ntheta,3))
  col = c()
  index = c()
  k = 0
  for(z in seq(0,1,len=6)) {
    index = c(index,k+(1:ntheta)-1,-1)
    for(theta in seq(0,1,len=ntheta)) {
      x = z*cos(2*pi*theta)/2 + 0.5
      y = z*sin(2*pi*theta)/2 + 0.5
      k = k + 1
      points[k,] = c(x,y,z)
      col[k] = hsv(theta,1,z)
    }
  }
  cat("Shape { appearance ",file=vrml.file);
  vrml.material.wrl(vrml.color(grey(0.5)),file=vrml.file)
  cat("geometry IndexedLineSet {",
      "coord Coordinate {\n",
      "point [",t(points),"]\n }\n",
      "colorPerVertex TRUE\n",
      "color Color { color [",t(vrml.color(col)),"]\n}",
      "coordIndex [",t(index),"]\n}}\n",file=vrml.file)
}


#' Depict colors geometrically
#'
#' Plots named colors as points in a three-dimensional cone.
#'
#' The colors are mapped into Hue-Saturation-Lightness (HSL).  Hue gives the
#' angle around the cone, Saturation the radial distance from the center line,
#' and Lightness the height from the base of the cone.
#'
#' @param col a vector of strings, naming colors.
#' @param cex a number controlling the size of the points.
#' @param light If \code{TRUE}, the cone will be surrounded by point light
#' sources.
#' @return Produces a VRML file which is opened by a VRML viewer.
#' @note The high-lightness, high-saturation part of the cone will always be
#' empty, because these colors are not achievable on computer displays (due to
#' the RGB representation).
#' @author Tom Minka
#' @references Rich Franzen's Wheel of Saturation, Intensity, and Hue.
#' \url{http://home.att.net/~rocq/SIHwheel.html}
#'
#' Charles Poynton's Color FAQ.
#' \url{http://www.poynton.com/notes/colour_and_gamma/ColorFAQ.html}
#' @examples
#'
#' color.cone(YlGnBu.colors(8))
#' color.cone(YR.colors(16))
#' color.cone(RYB.colors(7))
#' color.cone(topo.colors(20))
#' # reveals how topo.colors is not sequential in lightness
#'
color.cone <- function(col,cex=2,light=T,...) {
  # convert to HSY
  h = col2hsv(col)
  # replace Value by Lightness
  r = col2rgb(col)/255
  if(F) {
    r = r^3
    h[3,] = array(c(0.2126,0.7152,0.0722),c(1,3)) %*% r
    # lightness
    h[3,] = h[3,]^(1/3)
  } else {
    # approximation to lightness
    h[3,] = sqrt(colSums(r*r)/3)
  }
  h[3,] = h[3,]*255
  r <- hsv2cone(h)/255
  r = as.data.frame(t(r))
  x = r$x
  y = r$y
  z = r$z
  xlim = c(-1,1)
  ylim = c(-1,1)
  zlim = c(0,1)
  vrml.open(xlim,ylim,zlim,...)
  if(light) {
    cat("PointLight { ambientIntensity 0.5 intensity 0.5 location 0 0 0 }\n",
        file=vrml.file)
    cat("PointLight { intensity 0.5 location 1 1 0 }\n",
        file=vrml.file)
    cat("PointLight { intensity 0.5 location 0 1 0 }\n",
        file=vrml.file)
    cat("PointLight { intensity 0.5 location 1 0 0 }\n",
        file=vrml.file)
  }
  cat("Group { children [\n",file=vrml.file)
  vrml.axes(xlim,ylim,zlim,"x","y","z")
  if(T) {
    # outline of the cone
    cat("Transform { translation",cbind(0.5,0.5,0.5),"\n",
        "rotation",cbind(1,0,0,-pi/2),"\n",
        "scale",cbind(1,1,1)/2,"\n",
        "children Shape { ",file=vrml.file)
    cat("appearance ",file=vrml.file)
    vrml.material.wrl(diffuseColor=vrml.color(1),transparency=0.9,
                      specularColor=vrml.color(8),shininess=1,
                      file=vrml.file)
    cat("geometry Cone {}",file=vrml.file)
    cat("}\n}\n",file=vrml.file)

    # grey line
    vrml.line.wrl(cbind(0.5,0.5,0),cbind(0.5,0.5,1))
    vrml.cone.grid()
  }
  cat("]}\n",file=vrml.file)
  cat("Group { children [\n",file=vrml.file)
  vrml.sphere(x,y,z,col=col,cex=cex,...)
  p1 = vrml.to.scale(cbind(x[1],y[1],z[1]))
  n = length(x)
  for(i in 2:n) {
    p2 = vrml.to.scale(cbind(x[i],y[i],z[i]))
    vrml.line.wrl(p1,p2)
    p1 = p2
  }
  cat("]}\n",file=vrml.file)
  vrml.close()
}
test.color.cone <- function() {
  x = expand.grid(seq(0,1),seq(0,1),seq(0,1))
  col = rgb(x[,1],x[,2],x[,3])
  color.cone(col)
}
