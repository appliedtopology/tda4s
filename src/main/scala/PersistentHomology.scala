package org.appliedtopology.tda4s

import collection.mutable
import scala.collection.immutable.{SortedMap, SortedSet}
import scala.collection.Set

trait SimplexStream:
  def simplices(): Iterator[Simplex]
  def filtrationValues: PartialFunction[Simplex, Double]
  def simplicesInDimension(d: Int): Iterator[Simplex]

object PersistentHomology:
  def persistentHomology[T: Field as field](
      simplexStream: Iterator[Simplex],
      filtrationValues: PartialFunction[Simplex, Double]
  ): Seq[(Int, Double, Double)] =
    NaivePersistentHomology(simplexStream, filtrationValues)(using field).barcode
end PersistentHomology

class NaivePersistentHomology[T: Field as field](simplexstream: Iterator[Simplex], filtrationValues: PartialFunction[Simplex, Double]):
  import field._
  def filtrationValue(simplex: Simplex): Double =
    filtrationValues.applyOrElse(simplex, _ => Double.PositiveInfinity)
  def streamOrdering: Ordering[Simplex] = Ordering.by(filtrationValue).orElseBy(_.vertices)(Ordering.Implicits.sortedSetOrdering)
  case class State(
      cycleBasis: Map[Simplex, Chain[T]],
      boundaryBasis: Map[Simplex, Chain[T]],
      coboundaryLookup: Map[Simplex, Chain[T]],
      cycleBirths: Map[Simplex, Double]
  )
  var state: State = State(Map.empty, Map.empty, Map.empty, Map.empty)
  val barcodes: mutable.ListBuffer[(Int, Double, Double)] = mutable.ListBuffer.empty[(Int, Double, Double)]
  def step(simplex: Simplex): Unit = if (simplex.dimension <= 0) {
    val birth = filtrationValue(simplex)
    state = state.copy(
      cycleBasis = state.cycleBasis + (simplex -> Chain(
        SortedMap(simplex -> one)(using streamOrdering)
      )),
      cycleBirths = state.cycleBirths + (simplex -> birth)
    )
  } else {
    val simplexBoundary: Chain[T] = simplex.boundary(streamOrdering)(using field)
    val (boundaryModBoundaries, boundaryReductionLog) = Chain.reduce(simplexBoundary, state.boundaryBasis)
    if (boundaryModBoundaries.isZero) {
      // then ∂s is already a boundary, say the boundary of z
      // hence ∂(s-z) = ∂s - ∂z = 0
      // so s-z is a new cycle
      val cycle = // look through reduction reductionHistory, read off coboundaries and scalingFactor and add
        boundaryReductionLog
          .foldLeft(Chain.from(simplex)(streamOrdering)) { case (chain, (spx, coeff)) =>
            chain.scaleAdd(chain.field.neg(coeff.asInstanceOf[chain.field.F]), state.coboundaryLookup(spx))
          }

      val (cycleModCycles, _) = Chain.reduce(cycle, state.cycleBasis)
      val birthTime = filtrationValue(simplex)

      state = cycleModCycles.leading match {
        case Some((spx, _)) =>
          state.copy(
            cycleBasis = state.cycleBasis + (spx -> cycleModCycles),
            cycleBirths = state.cycleBirths + (spx -> birthTime)
          )
        case None: Option[(Simplex, Double)] => state
      }
    } else {
      // since ∂s is not a boundary already, it has to be a cycle
      // reducing mod cycles will tell us which cycle just became a boundary with the addition of simplex
      // val (dsModZ, reductionZ) = Chain.reduceBy(boundaryModBoundaries, state.cycleBasis)
      // assert(dsModZ.isZero)
      val dyingCycleLeadingSimplex = boundaryModBoundaries.simplexCoefficients.keys
        .filter(k => state.cycleBasis.contains(k))
        .max(using streamOrdering)
      val updatedCoboundary = // look through reduction reductionHistory, read off coboundaries and scalingFactor and add
        boundaryReductionLog
          .foldLeft(Chain.from(simplex)(streamOrdering)) { case (chain, (spx, coeff)) =>
            chain.scaleAdd(chain.field.neg(coeff.asInstanceOf[chain.field.F]), state.coboundaryLookup(spx))
          }
      if (math.abs(filtrationValue(simplex) - state.cycleBirths(dyingCycleLeadingSimplex)) > 1e-15)
        barcodes.append((dyingCycleLeadingSimplex.dimension, state.cycleBirths(dyingCycleLeadingSimplex), filtrationValue(simplex)))
      assert(simplex.dimension == dyingCycleLeadingSimplex.dimension + 1)
      state = state.copy(
        cycleBasis = state.cycleBasis - dyingCycleLeadingSimplex,
        boundaryBasis = state.boundaryBasis + (boundaryModBoundaries.leading.get._1 -> boundaryModBoundaries),
        cycleBirths = state.cycleBirths - dyingCycleLeadingSimplex,
        coboundaryLookup = state.coboundaryLookup + (dyingCycleLeadingSimplex -> updatedCoboundary)
      )
    }
  }
  def barcode: Seq[(Int, Double, Double)] = {
    simplexstream.foreach(step)
    state.cycleBirths
      .map { case (simplex, birth) =>
        (simplex.dimension, birth, Double.PositiveInfinity)
      }
      .toSeq
      .concat(barcodes)
  }
end NaivePersistentHomology
