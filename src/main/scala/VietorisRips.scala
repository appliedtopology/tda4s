package org.appliedtopology.tda4s

import scala.collection.Set

/** Vietoris Rips trait that work as an interface for generate Vietoris Rips simplex streams from classes that inherit it. It provides a
  * default implementation for Vietoris-Rips filtration values
  *
  * @param metricSpace
  *   metric space where the underlying metric is used to compute distances between points for the filtration values
  */
trait VietorisRips(val metricSpace: MetricSpace[Int]) extends SimplexStream:
  // partial function that computes the filtration of a simplex
  // returns -infinity is the simplex is empty,
  // 0 for a 0-dimensional simplex, and the largest pairwise distance for a higher dimensional simplex
  def filtrationValues: PartialFunction[Simplex, Double] = {
    case spx if spx.vertices.isEmpty => Double.NegativeInfinity
    case spx if spx.dimension == 0   => 0.0
    case spx =>
      spx.vertices.toSeq.combinations(2).map { case Seq(x, y) => metricSpace.distance(x, y) }.max
  }

  // Iterator for simplexes
  def simplicesInDimension(dimension: Int): Iterator[Simplex]

  // Get the simplexes by iterating through the elements of the metric space and storing them in a flattened list
  def simplices(): Iterator[Simplex] = (0 to metricSpace.elements.size).iterator.flatMap(simplicesInDimension)

/** Class for the Naive Vietoris Rips Stream formulation. Naive in the sense that no optimizations are applied to speed computation. It
  * extends the Vietoris Rips Trait
  *
  * @param metricSpace
  *   metric space where the underlying metric is used to compute distances between points for the filtration values
  */
class NaiveVietorisRips(metricSpace: MetricSpace[Int]) extends VietorisRips(metricSpace):
  /** Iterator object to generate all of the simplexes of a given dimension. Generates all of the combinations of a given dimension and
    * computes the resulting filtration value by taking the largest pairwise distance and sorting accordingly
    *
    * @param metricSpace
    *   metric space where the underlying metric is used to compute distances between points for the filtration values
    */
  override def simplicesInDimension(dimension: Int): Iterator[Simplex] =
    metricSpace.elements
      .combinations(dimension + 1)
      .map(Simplex.from(_))
      .toSeq
      .sortBy(filtrationValues.applyOrElse(_, _ => Double.PositiveInfinity))
      .iterator
