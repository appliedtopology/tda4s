package org.appliedtopology.tda4s

import collection.mutable
import scala.collection.immutable.{SortedMap, SortedSet}
import scala.collection.Set

trait SimplexStream:
  def simplices(): Iterator[Simplex]
  def filtrationValues: PartialFunction[Simplex, Double]
  def simplicesInDimension(d: Int): Iterator[Simplex]

object PersistentHomology:
  // Companion object for the NaivePersistentHomology class
  def persistentHomology[T: Field as field](
      simplexStream: Iterator[Simplex],
      filtrationValues: PartialFunction[Simplex, Double]
  ): Seq[(Int, Double, Double)] =
    NaivePersistentHomology(simplexStream, filtrationValues)(using field).barcode
end PersistentHomology

/** Class to run the persistent homology algorithm. It is naive in the sense that it is a textbook implementation with no optimizations
  * applied to it. It utilizes a form of row reduction to find the pivots of the boundary matrix, which correspond to cycles in the
  * underlying topology.
  *
  * @tparam Field[T]
  *   type parameter for the field, requires we have an instantiation of F[T]
  *
  * @param simplexStream
  *   iterable object of simplexes
  *
  * @param filtrationValues
  *   partial function that takes a simplex as an input and returns the filtration value as a double
  */
class NaivePersistentHomology[T: Field as field](simplexstream: Iterator[Simplex], filtrationValues: PartialFunction[Simplex, Double]):
  import field._

  /** Method for the filtration value of a simplex. Returns either the filtation value specified in filtrationValues or infinity if it is
    * not defined
    *
    * @param simplex
    *   The simplex used to compute the filtration value
    */
  def filtrationValue(simplex: Simplex): Double =
    filtrationValues.applyOrElse(simplex, _ => Double.PositiveInfinity)

  /** Method to order the simplexes in the simplex stream. Should occur by the filtration values or the default set ordering in Scala if the
    * filtration isn't defined.
    */
  def streamOrdering: Ordering[Simplex] = Ordering.by(filtrationValue).orElseBy(_.vertices)(Ordering.Implicits.sortedSetOrdering)

  // Homology state is tracked through a case class with lookup tables for the cycle basis, boundary basis, coboundary, and cycle births
  case class State(
      cycleBasis: Map[Simplex, Chain[T]],
      boundaryBasis: Map[Simplex, Chain[T]],
      coboundaryLookup: Map[Simplex, Chain[T]],
      cycleBirths: Map[Simplex, Double]
  )

  // Homology state is null by default
  var state: State = State(Map.empty, Map.empty, Map.empty, Map.empty)
  // barcodes are implemented with a mutable (linked) list and each entry is (Dimension, birth, death)
  val barcodes: mutable.ListBuffer[(Int, Double, Double)] = mutable.ListBuffer.empty[(Int, Double, Double)]

  /** Method to perform a step in the persistent homology computation using a simplex from a stream
    *
    * @param simplex
    *   Current simplex to run the computation
    */
  def step(simplex: Simplex): Unit =
    if (simplex.dimension <= 0) {
      val birth = filtrationValue(simplex)
      state = state.copy(
        cycleBasis = state.cycleBasis + (simplex -> Chain(
          SortedMap(simplex -> one)(using streamOrdering)
        )),
        cycleBirths = state.cycleBirths + (simplex -> birth)
      )
      // Degree 0 simplexes are straightforward for the computation and everything after the next else statement is for higher-dimensional simplexes
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
  // updates the homology state and generates barcode by applying the step function
  // to each simplex in the simplex stream
end NaivePersistentHomology
