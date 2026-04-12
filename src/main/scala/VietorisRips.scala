package org.appliedtopology.tda4s

import scala.collection.Set

trait VietorisRips(val metricSpace: MetricSpace[Int]) extends SimplexStream:
  def filtrationValues: PartialFunction[Simplex, Double] = {
    case spx if spx.vertices.isEmpty => Double.NegativeInfinity
    case spx if spx.dimension == 0   => 0.0
    case spx =>
      spx.vertices.toSeq.combinations(2).map { case Seq(x, y) => metricSpace.distance(x, y) }.max
  }
  def simplicesInDimension(dimension: Int): Iterator[Simplex]
  def simplices(): Iterator[Simplex] = (0 to metricSpace.elements.size).iterator.flatMap(simplicesInDimension)

class NaiveVietorisRips(metricSpace: MetricSpace[Int]) extends VietorisRips(metricSpace):
  override def simplicesInDimension(dimension: Int): Iterator[Simplex] =
    metricSpace.elements
      .combinations(dimension + 1)
      .map(Simplex.from(_))
      .toSeq
      .sortBy(filtrationValues.applyOrElse(_, _ => Double.PositiveInfinity))
      .iterator
