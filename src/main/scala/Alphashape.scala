package org.appliedtopology.tda4s

import collection.mutable
import scala.math.{pow, sqrt}
import org.apache.commons.math3.linear.{MatrixUtils, QRDecomposition, RealMatrix, RealVector, SingularValueDecomposition}
import com.gurobi.gurobi.*
import io.github.zabuzard.closy.external.*
import org.apache.commons.math3.linear.MatrixUtils.createRealVector

import scala.jdk.CollectionConverters.*

def Alpha(pts: Seq[Array[Double]]): SimplexStream = pts match {
  case pts if pts.isEmpty => HelixDelaunay(pts)
  case pts if pts.head.size > 0 =>
    HelixDelaunay(pts) // Helix should be faster for dim: 7 - 17. Adjust this check when additional impl exists.
}


// utilities for Delaunay computations
type Point = RealVector
object Point:
  def apply(coords : Array[Double]): Point = MatrixUtils.createRealVector(coords)

case class Hyperplane(normal: Point, offset: Double):
  def isLight(point: Point): Boolean = normal.dotProduct(point) >= offset

  def dist(point: Point): Double = normal.dotProduct(point) - offset

  def reverse: Hyperplane = Hyperplane(normal.mapMultiply(-1.0), -offset)

object Hyperplane:
  def from(pointSet: Seq[Point]): Hyperplane = {
    assert(pointSet.size >= pointSet.head.getDimension)
    val points = pointSet.map(p => p.subtract(pointSet.head))
    // now we need to find a basis of orthogonal complement to the space spanned by these
    val A = MatrixUtils.createRealMatrix(points.toSeq.map(_.toArray).toArray)
    val svdA = SingularValueDecomposition(A.transpose)
    val nullSV = svdA.getSingularValues.zipWithIndex.filter(_._1.abs < 1e-5)
    val n = nullSV.map((s, i) => svdA.getU.getColumnVector(i)).head // should leave out head if we want higher codimensions
    val offset = n.dotProduct(pointSet.head)
    Hyperplane(n, offset)
  }

case class Ridge(normals : (Point,Point))
object Ridge:
  def apply(pointSet : Seq[Point], hyperplane: Hyperplane): Ridge = {
    assert(pointSet.size == hyperplane.normal.getDimension-1)
    val Hyperplane(normal1, offset1) = hyperplane
    val points = pointSet.map(p => p.subtract(pointSet.head)).appended(normal1)
    val Hyperplane(normal2, offset) = Hyperplane.from(points)
    Ridge((normal1, normal2))
  }

case class Hypersphere(center: Point, radius: Double):
  def contains(point: Point): Boolean = point.subtract(center).getNorm < radius - 1e-5

  def onSphere(point: Point): Boolean = math.abs(point.subtract(center).getNorm - radius) < 1e-5

object Hypersphere:
  def apply(pointSet: Seq[Point]): Hypersphere = {
    val vvs = for {
      v1 <- pointSet.indices
      v2 <- pointSet.indices
      if v1 < v2
    } yield (pointSet(v1).subtract(pointSet(v2)), math.pow(pointSet(v1).getNorm, 2) - math.pow(pointSet(v2).getNorm, 2))
    val V = MatrixUtils.createRealMatrix(vvs.map(_._1.toArray).toArray)
    val b = MatrixUtils.createRealVector(vvs.map(_._2).toArray).mapMultiply(0.5)
    val c = SingularValueDecomposition(V).getSolver.solve(b)
    new Hypersphere(c, c.getDistance(pointSet.head))
  }

case class DelaunaySimplex(simplex: Simplex, circumsphere: Hypersphere)

/**
 * Giftwrapping algorithm in n+1 dimensions to find Delaunay triangulation with lower convex hull
 */
class GiftwrappingDelaunay(pts: Seq[Array[Double]]) {
  val points : Seq[Point] = pts.map(Point.apply)
  val copoints = points.map(p => p.append(p.getNorm*p.getNorm))
  val ambientDimension = points.head.getDimension

  // find seed facet on convex hull
  var seedFacet = points.indices.sortBy(i => points(i).getEntry(0)).take(ambientDimension+1).toSeq
  var seedHP = Hyperplane.from(seedFacet.map(copoints))

}

/** Based on https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=10917453&tag=1
  * @param pts
  */
class HelixDelaunay(pts: Seq[Array[Double]]) extends SimplexStream {
  val points: Seq[Point] = pts.map(MatrixUtils.createRealVector(_))
  val ambientDimension: Int = points.head.getDimension

  val validated: mutable.Set[DelaunaySimplex] = mutable.Set.empty
  val candidate: mutable.ArrayDeque[Simplex] = mutable.ArrayDeque.empty

  var startingSimplex: Set[Int] = points.indices.take(ambientDimension).toSet
  {
    var convexHullHP: Boolean = false
    while (!convexHullHP) {
      val currentHP: Hyperplane = Hyperplane.from(startingSimplex.map(p => points(p)).toSeq)
      val lightPoints = (points.indices.toSet -- startingSimplex).filter(pi => currentHP.dist(points(pi)) > 1e-5)
      if (lightPoints.isEmpty) convexHullHP = true
      else {
        startingSimplex = (lightPoints ++ startingSimplex).toSeq
          .sortBy(_ => util.Random.nextDouble())
          .take(ambientDimension)
          .toSet
        //      startingSimplex = (lightPoints ++ startingSimplex)
        //        .map(p => p -> currentHP.dist(points(p)))
        //        .toSeq
        //        .sortBy(_._2)
        //        .takeRight(ambientDimension)
        //        .map(_._1)
        //        .toSet
      }
    }
    startingSimplex = points.indices.filter(pi => Hyperplane.from(startingSimplex.map(points).toSeq).dist(points(pi)).abs < 1e-5) match {
      case vs if vs.size == ambientDimension => vs.toSet
      case vs if vs.size > ambientDimension => Set(vs.head, vs.tail.minBy(vi => points(vs.head).getDistance(points(vi))))
    }

    // brute force search for first delaunay simplex
    var done = false
    for (pi <- points.indices)
      if (!done) {
        var valid = true
        if (!startingSimplex.contains(pi)) {
          val circumsphere = Hypersphere((startingSimplex + pi).map(p => points(p)).toSeq)
          val containedPoints = points.indices.toSet.filter(qi => circumsphere.contains(points(qi)))
          val containedDistances = containedPoints.map(qi => circumsphere.center.getDistance(points(qi)))
          if (
            (points.indices.toSet -- (startingSimplex + pi))
              .exists(qi => circumsphere.contains(points(qi)) && !circumsphere.onSphere(points(qi)))
          ) {
            valid = false
          }
          if (valid) {
            done = true
            validated.add(DelaunaySimplex(Simplex.from((startingSimplex + pi).toSeq), circumsphere))
          }
        }
      }
    assert(validated.nonEmpty)
  }
  case class FrontierCase(facet: Simplex, hyperplane: Hyperplane, complement: Int, cofacetHypersphere: Hypersphere)

  object FrontierCase:
    def apply(cofacet: DelaunaySimplex, complement: Int): FrontierCase = {
      val facet = cofacet.simplex - complement
      val hyperplane = Hyperplane.from(facet.vertices.toSeq.map(points(_)))
      if (hyperplane.isLight(points(complement)))
        FrontierCase(facet, hyperplane.reverse, complement, cofacet.circumsphere)
      else
        FrontierCase(facet, hyperplane, complement, cofacet.circumsphere)
    }

  val frontierCases: mutable.ArrayDeque[FrontierCase] = mutable.ArrayDeque.empty
  val visitedFacets: mutable.Set[Simplex] = mutable.Set(Simplex.from(startingSimplex))
  val cospherical: mutable.Set[Set[Int]] = mutable.Set.empty
  println(s"starting facet: $startingSimplex")

  def addFrontierCase(simplex: DelaunaySimplex, complement: Int): Unit =
    frontierCases.zipWithIndex.find((fc: FrontierCase, i: Int) => fc.facet.vertices == simplex.simplex.vertices) match {
      case None =>
        frontierCases.addOne(
          FrontierCase(simplex, complement)
        )
      case Some((otherCase, otherIndex)) => frontierCases.remove(otherIndex)
    }

  def handleCosphericalPoints(cosphericalPoints: Seq[Int], frontierCase: FrontierCase, newDelaunaySimplex: DelaunaySimplex): Unit = {
    val spherepoints: mutable.SortedSet[Int] = cosphericalPoints.to(mutable.SortedSet)
    // all of these work as extra point: every subset of the spherepoints is a valid Delaunay simplex with this empty circumsphere
    // so we need to pick a tiling subset of them. We start a local version of this frontier walking algorithm
    // we also need to make sure we don't come back inside this cospherical point set in a later iteration
    cospherical.add(spherepoints.toSet)
    println(s"o   cospherical: $cospherical")
    frontierCases.removeAll((fc: FrontierCase) => fc.facet.vertices.subsetOf(frontierCase.facet.vertices + frontierCase.complement))
    spherepoints.subtractAll(newDelaunaySimplex.simplex.vertices)
    val facets: mutable.ArrayDeque[(Simplex, Simplex)] =
      mutable.ArrayDeque.from(frontierCase.facet.vertices.toSeq.map(fi => (newDelaunaySimplex.simplex - fi, newDelaunaySimplex.simplex)))
    while (facets.nonEmpty) {
      // take a facet
      val (facet, cofacet) = facets.removeHead()
      val hyperplane: Hyperplane = Hyperplane.from(facet.vertices.toSeq.map(points)) match {
        case candidateHP if frontierCase.facet.vertices.toSeq.map(points).exists(candidateHP.isLight) => candidateHP.reverse
        case candidateHP                                                                              => candidateHP
      } // TODO: this doesn't work, maybe the Deque is the wrong data structure here?
      spherepoints.filter(hyperplane.isLight.compose(points)).headOption match {
        case Some(pi) =>
          val nds = DelaunaySimplex(facet + pi, newDelaunaySimplex.circumsphere)
          validated.add(nds)
          println(s"o +   added $nds")
          facets.addAll(facet.vertices.toSeq.map(pj => (Simplex(nds.simplex.vertices - pj), nds.simplex)))
          spherepoints.remove(pi)
        case None =>
          addFrontierCase(DelaunaySimplex(cofacet, newDelaunaySimplex.circumsphere), cofacet.vertices.diff(facet.vertices).head)
      }
    }
  }

  val seedDelaunaySimplex = validated.head
  points.indices
    .map(i => (i, seedDelaunaySimplex.circumsphere.center.getDistance(points(i))))
    .filter((i, d) => math.abs(d - seedDelaunaySimplex.circumsphere.radius) <= 1e-5)
    .map(_._1)
    .to(mutable.SortedSet) match {
    case spherepoints if spherepoints.size > ambientDimension + 1 =>
      handleCosphericalPoints(
        spherepoints.toSeq,
        FrontierCase(seedDelaunaySimplex, (seedDelaunaySimplex.simplex.vertices -- startingSimplex).head),
        seedDelaunaySimplex
      )
    case spherepoints => frontierCases.addAll(startingSimplex.map(pi => FrontierCase(seedDelaunaySimplex, pi)).toSeq)
  }

  // handle a frontier case
  while (frontierCases.nonEmpty) {
    val frontierCase = frontierCases.removeHead()
    if (visitedFacets.contains(frontierCase.facet)) {
      println(s"x   skipping $frontierCase")
    } else {
      println(frontierCase)
      visitedFacets.add(frontierCase.facet)
      // "light points" (in front of the facet) split into inside and outside a small circumsphere of the facet
//      val circumsphere = Hypersphere(frontierCase.facet.vertices.toSeq.map(points))
//      val (inside, outside) = (points.indices.toSet -- frontierCase.facet.vertices)
//        .filter(pi => frontierCase.hyperplane.isLight(points(pi)))
//        .partition(pi => circumsphere.contains(points(pi)))
//      if ((inside ++ outside).nonEmpty) {
//        val newDelaunaySimplex: DelaunaySimplex = if (inside.nonEmpty) {
//          inside
//            .map(pi => frontierCase.facet.vertices + pi)
//            .map(vs => DelaunaySimplex(Simplex(vs), Hypersphere(vs.toSeq.map(points))))
//            .maxBy(_.circumsphere.radius)
//        } else {
//          outside
//            // .filter(pi => frontierCase.cofacetHypersphere.center.getDistance(points(pi)) <= 2 * frontierCase.cofacetHypersphere.radius)
//            .map(pi => frontierCase.facet.vertices + pi)
//            .map(vs => DelaunaySimplex(Simplex(vs), Hypersphere(vs.toSeq.map(points))))
//            .minBy(_.circumsphere.radius)
//        }
      (points.indices.toSet -- frontierCase.facet.vertices).toSeq
        .filter(pi => frontierCase.hyperplane.isLight(points(pi)))
        .sortBy(pi => frontierCase.cofacetHypersphere.center.getDistance(points(pi)))
        .view
        .map ( pi => DelaunaySimplex(frontierCase.facet + pi, Hypersphere((frontierCase.facet + pi).vertices.toSeq.map(points))) )
        .collectFirst { case ds if (!points.exists(ds.circumsphere.contains)) => ds } match {
        case Some(newDelaunaySimplex) : Option[DelaunaySimplex] if !validated.exists(ds => newDelaunaySimplex.simplex == ds.simplex) =>
          // check whether we have "too many" cospherical points; in that case we have to tile them on our own
          val spherepoints: mutable.SortedSet[Int] = points.indices
            .map(i => (i, newDelaunaySimplex.circumsphere.center.getDistance(points(i))))
            .filter((i, d) => math.abs(d - newDelaunaySimplex.circumsphere.radius) <= 1e-5)
            .map(_._1)
            .to(mutable.SortedSet)
          if ((spherepoints.size > ambientDimension + 1)) {
            if (!cospherical.contains(spherepoints.toSet)) {
              handleCosphericalPoints(spherepoints.toSeq, frontierCase, newDelaunaySimplex)
            }
          } else {
            validated.add(newDelaunaySimplex)
            println(s"  +   added $newDelaunaySimplex")
            frontierCase.facet.vertices.toSeq
              .foreach(vi => addFrontierCase(newDelaunaySimplex, vi))
          }
        case Some(newDelaunaySimplex) : Option[DelaunaySimplex] => ()
        case None => ()
      }
    }
  }

  val simplicesMap: Map[Int, Seq[Simplex]] = Map.from(
    (0 to ambientDimension).map(d => d -> validated.map(ds => ds.simplex.vertices.subsets(d + 1)).flatMap(_.toSet).map(Simplex.apply).toSeq)
  )

  def edgeIsDelaunay(s: Simplex): Option[Double] = {
    val Seq(src, tgt) = s.vertices.toSeq
    val circumcenter: Point = points(src).add(points(tgt)).mapMultiply(0.5)
    val circumsphere = Hypersphere(circumcenter, circumcenter.getDistance(points(src)))
    if (points.filter(p => circumcenter.getDistance(p) < circumsphere.radius - 1e-5).size > 0) None
    else Some(circumsphere.radius)
  }

  def computeFVal(s: Simplex): Double = s.dimension match {
    case 0 => 0.0
    case 1 =>
      edgeIsDelaunay(s) match {
        case Some(alpha) => alpha
        case None =>
          validated
            .filter(v => s.vertices.subsetOf(v.simplex.vertices))
            .map(ds => ds.circumsphere.radius)
            .min
      }
    case _ =>
      validated
        .filter(v => s.vertices.subsetOf(v.simplex.vertices))
        .map(ds => ds.circumsphere.radius)
        .min
  }

  override def filtrationValues: PartialFunction[Simplex, Double] =
    Map.from[Simplex, Double](
      (0 to ambientDimension)
        .flatMap(simplicesMap)
        .map(s => s -> computeFVal(s))
    )

  val simplicesSortedMap: Map[Int, Seq[Simplex]] = simplicesMap.map((d, v) => (d, v.sortBy(filtrationValues)))

  override def simplicesInDimension(d: Int): Iterator[Simplex] = simplicesSortedMap(d).iterator

  override def simplices(): Iterator[Simplex] = (0 to ambientDimension).iterator.flatMap(simplicesInDimension)
}

class QuadProgAlpha(pts: Seq[Array[Double]], maxDistance: Double = Double.PositiveInfinity) extends SimplexStream:
  type Point = Array[Double]
  val points: Seq[Point] = pts

  class EuclideanDistance extends Metric[Point] {
    override def distance(first: Point, second: Point): Double =
      sqrt(first.zip(second).map((x0, x1) => pow(x0 - x1, 2)).sum)
  }
  val pointMetric = EuclideanDistance()
  class EuclideanPointDistance extends Metric[Int] {
    override def distance(first: Int, second: Int): Double =
      pointMetric.distance(points(first), points(second))
  }
  val metric = EuclideanPointDistance()
  var nnc = NearestNeighborComputations.of(metric)
  points.indices.foreach(nnc.add(_))

  val maxd: Double = if (maxDistance == Double.PositiveInfinity) {
    points.indices.map(p => points.indices.map(q => metric.distance(p, q)).max).min
  } else maxDistance

  val filtrationValues: mutable.Map[Simplex, Double] = mutable.Map(points.indices.map(Simplex(_) -> 0.0)*)
  val simplexMap: mutable.Map[Int, Seq[Simplex]] = mutable.Map(0 -> points.indices.map(Simplex(_)))
  val graph: Map[Int, mutable.Set[Int]] = points.indices.map(i => i -> mutable.Set.empty).toMap
  override def simplices(): Iterator[Simplex] = (0 to points.size).iterator.flatMap(simplicesInDimension)
  def simplicesInDimension(dimension: Int): Iterator[Simplex] = dimension match {
    case i if simplexMap.contains(i) => simplexMap(i).iterator
    case i if i > points.head.size   => Iterator.empty
    case 1 =>
      val output: mutable.Set[Simplex] = mutable.Set.empty
      // first off, nearest neighbor is always a Delaunay edge
      for (x <- points.indices) {
        val nn = nnc.getKNearestNeighbors(x, 2).asScala.last
        val simplex = Simplex.from(Seq(x, nn))
        graph(nn).add(x)
        graph(x).add(nn)
        output.add(simplex)
        filtrationValues(simplex) = metric.distance(x, nn) / 2
      }
      // next up, if the circle with diameter given by the pair of points is empty, it's a Delaunay edge
      // because the Gabriel graph includes into the Delaunay graph
      for (x <- points.indices) {
        val env = GRBEnv(true)
        env.set(GRB.IntParam.OutputFlag, 0)
        env.start()
        val model = GRBModel(env)

        val ys = model.addVars(points(x).size, GRB.CONTINUOUS)

        val obj = GRBQuadExpr()
        // |x-y|^2 = (x0-y0)^2 + (x1-y1)^2 + ... + (xd-yd)^2
        //    = x0^2 - 2x0y0 + y0^2 + ....
        //    = sum xi^2 -2 sum xiyi + sum yi^2
        // is maximized when the variable bits are maximized
        //    ~ -2 sum xi yi + sum yi^2
        ys.indices.foreach { i =>
          obj.addConstant(pow(points(x)(i), 2))
          obj.addTerm(1.0, ys(i), ys(i))
          obj.addTerm(-2.0 * points(x)(i), ys(i))
        }
        model.setObjective(obj, GRB.MINIMIZE)

        // for each other point, we have similar-ish requirements:
        // |x-y|^2 = or ≤ |z-y|^2
        // sum xi^2 - 2 sum xi yi + sum yi^2 =/≥ sum zi^2 - 2sum zi yi + sum yi^2
        // cancel the yi^2
        // sum xi^2 - sum zi^2 - 2 sum yi (xi-zi)
        // these are EQUAL for vertices in a Delaunay/Alpha simplex, and ≤ for non-simplex vertices.
        points.indices.filter(_ != x).foreach { z =>
          val constr = GRBLinExpr() // linear in the yi
          ys.indices.foreach { i =>
            constr.addConstant(pow(points(x)(i), 2) - pow(points(z)(i), 2))
            constr.addTerm(-2.0 * (points(x)(i) - points(z)(i)), ys(i))
          }
          model.addConstr(constr, GRB.LESS_EQUAL, 0.0, z.toString)
        }
        model.update()

        for (y <- points.indices)
          if ((x < y) && !graph(x).contains(y)) {

            // the Gabriel criterion: the hypersphere with radius d around m must not contain any other points
            val m: Point = points(x).zip(points(y)).map((xi, yi) => (xi + yi) / 2)
            val d = pointMetric.distance(m, points(x))
            if (
              points.indices
                .filter(i => (i != x) && (i != y))
                .map(points(_))
                .forall(p => pointMetric.distance(m, p) >= d - 1e-5)
            ) {
              val simplex = Simplex.from(Seq(x, y))
              graph(y).add(x)
              graph(x).add(y)
              output.add(simplex)
              filtrationValues(simplex) = d
            } else {
              val yconstr = model.getConstrByName(y.toString)
              yconstr.set(GRB.CharAttr.Sense, GRB.EQUAL)
              model.optimize()
              println(s"[$x,$y]: Optimizer status ${model.get(GRB.IntAttr.Status)}")
              if (model.get(GRB.IntAttr.Status) == GRB.OPTIMAL) {
                // solution found, so this is a Delaunay edge
                val simplex = Simplex.from(Seq(x, y))
                graph(y).add(x)
                graph(x).add(y)
                output.add(simplex)
                filtrationValues(simplex) = d
              }
              yconstr.set(GRB.CharAttr.Sense, GRB.LESS_EQUAL)
            }
          }

        model.dispose()
        env.dispose()
      }
      simplexMap(1) = output.toSeq.sortBy(filtrationValues)
      simplexMap(1).iterator
    case _ =>
      val facets = simplicesInDimension(dimension - 1).toSeq
      val output: mutable.Set[Simplex] = mutable.Set.empty
      for (x <- points.indices) {
        val candidateFacets = facets
          .filter(f => x < f.vertices.head)
          .filter(f => f.vertices.forall(n => (dimension <= 1) || graph(n).contains(x)))
        if (candidateFacets.nonEmpty) {
          val neighborhood = candidateFacets.flatMap(_.vertices).toSet

          val env = GRBEnv(true)
          env.set(GRB.IntParam.OutputFlag, 0)
          env.start()
          val model = GRBModel(env)
          // val c1 = maxd
          // model.set(GRB.DoubleParam.Cutoff, c1)
          model.set(GRB.IntParam.DualReductions, 0)

          val ys = model.addVars(points(x).size, GRB.CONTINUOUS)

          val obj = GRBQuadExpr()
          // |x-y|^2 = (x0-y0)^2 + (x1-y1)^2 + ... + (xd-yd)^2
          //    = x0^2 - 2x0y0 + y0^2 + ....
          //    = sum xi^2 -2 sum xiyi + sum yi^2
          // is maximized when the variable bits are maximized
          //    ~ -2 sum xi yi + sum yi^2
          ys.indices.foreach { i =>
            obj.addConstant(pow(points(x)(i), 2))
            obj.addTerm(1.0, ys(i), ys(i))
            obj.addTerm(-2.0 * points(x)(i), ys(i))
          }
          model.setObjective(obj, GRB.MINIMIZE)

          // for each other point, we have similar-ish requirements:
          // |x-y|^2 = or ≤ |z-y|^2
          // sum xi^2 - 2 sum xi yi + sum yi^2 =/≤ sum zi^2 - 2sum zi yi + sum yi^2
          // cancel the yi^2, move everything to the left
          // sum xi^2 - sum zi^2 - 2 sum yi (xi-zi) =/≤ 0
          // these are EQUAL for vertices in a Delaunay/Alpha simplex, and ≤ for non-simplex vertices.
          points.indices.filter(_ != x).foreach { z =>
            val constr = GRBLinExpr() // linear in the yi
            ys.indices.foreach { i =>
              constr.addConstant(pow(points(x)(i), 2) - pow(points(z)(i), 2))
              constr.addTerm(-2.0 * (points(x)(i) - points(z)(i)), ys(i))
            }
            model.addConstr(constr, GRB.LESS_EQUAL, 0.0, z.toString)
          }
          model.update()

          for (candidateFacet: Simplex <- candidateFacets) {
            for (z <- candidateFacet.vertices)
              model.getConstrByName(z.toString).set(GRB.CharAttr.Sense, GRB.EQUAL)
            model.update()

            model.optimize()
            if (Seq(GRB.OPTIMAL, GRB.LOCALLY_OPTIMAL, GRB.USER_OBJ_LIMIT).contains(model.get(GRB.IntAttr.Status))) {
              val fval = model.get(GRB.DoubleAttr.ObjVal)
              if (fval <= Double.PositiveInfinity) {
                filtrationValues(candidateFacet + x) = sqrt(2 * fval)
                output.add(candidateFacet + x)
                if (dimension == 1) {
                  graph(x).add(candidateFacet.vertices.head)
                  graph(candidateFacet.vertices.head).add(x)
                }
              }
            }
            for (z <- candidateFacet.vertices)
              model.getConstrByName(z.toString).set(GRB.CharAttr.Sense, GRB.LESS_EQUAL)
            model.update()
          }
          model.dispose()
          env.dispose()
        }
      }
      simplexMap(dimension) = output.toSeq.sortBy(filtrationValues)
      simplexMap(dimension).iterator
  }
