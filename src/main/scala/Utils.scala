package org.appliedtopology.tda4s

import scala.collection.immutable.IndexedSeqDefaults.defaultApplyPreferredMaxLength.{+, -}

object SVG:
  def plotComplex(simplices: Seq[Simplex], positions: Seq[Array[Double]]): String = {
    val xmin = positions.map(_(0)).min
    val ymin = positions.map(_(1)).min
    val xmax = positions.map(_(0)).max
    val ymax = positions.map(_(1)).max
    val width = xmax - xmin
    val height = ymax - ymin
    val viewPointXmin = xmin - 0.05*width
    val viewPointYmin = ymin - 0.05*height
    val viewPointWidth = 1.1*width
    val viewPointHeight = 1.1*height
    
    // Generate SVG elements for the simplices
    val svgElements = simplices.flatMap {
      case simplex if simplex.vertices.size == 1 =>
        // Single vertex -> add a small circle
        simplex.vertices.headOption.map { v =>
          val pos = positions(v)
          s"""<circle cx="${pos(0)}" cy="${pos(1)}" r="3" fill="black"><title>${v}</title></circle>"""
        }

      case simplex if simplex.vertices.size == 2 =>
        // Edge (line) -> add a line element
        for
          a <- simplex.vertices.headOption
          b <- simplex.vertices.toIndexedSeq.lift(1)
          pos1 = positions(a)
          pos2 = positions(b)
        yield
          s"""<line x1="${pos1(0)}" y1="${pos1(1)}" x2="${pos2(0)}" y2="${pos2(1)}" stroke="black" stroke-width="1"/>"""

      case simplex if simplex.vertices.size == 3 =>
        // Triangle -> add a polygon element
        for
          a <- simplex.vertices.headOption
          b <- simplex.vertices.toIndexedSeq.lift(1)
          c <- simplex.vertices.toIndexedSeq.lift(2)
        yield
          val pos1 = positions(a)
          val pos2 = positions(b)
          val pos3 = positions(c)
          s"""<polygon points="${pos1(0)},${pos1(1)} ${pos2(0)},${pos2(1)} ${pos3(0)},${pos3(1)}" fill="blue" fill-opacity="0.25" stroke="none"/>"""

      case _ =>
        // Ignore higher-dimensional simplices
        None
    }.mkString("\n")

    // Generate the full SVG element
    val svgTemplate =
      s"""
         |<svg xmlns="http://www.w3.org/2000/svg" viewBox="$viewPointXmin $viewPointYmin $viewPointWidth $viewPointHeight" width="500" height="500">
         |  $svgElements
         |</svg>
      """.stripMargin

    svgTemplate
  }
