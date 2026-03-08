package org.appliedtopology.tda4s

import collection.mutable
import scala.collection.immutable.SortedMap

object PersistentHomology:
  def persistentHomology(simplexstream: Seq[Set[Int]], filtrationValues: Map[Set[Int], Double]): Seq[(Int, Double, Double)] =
    NaivePersistentHomology(simplexstream, filtrationValues).barcode
end PersistentHomology

class NaivePersistentHomology(simplexstream: Seq[Set[Int]], filtrationValues: Map[Set[Int], Double]):
  case class State(
      cycleBasis: Map[Set[Int], Chain],
      boundaryBasis: Map[Set[Int], Chain],
      coboundaryLookup: Map[Set[Int], Set[Int]],
      cycleBirths: Map[Set[Int], Double]
  )
  var state = State(Map.empty, Map.empty, Map.empty, Map.empty)
  val barcodes = mutable.ListBuffer.empty[(Int, Double, Double)]
  def step(simplex: Set[Int]): Unit = if (simplex.size <= 1) {
    val birth = filtrationValues.getOrElse(simplex, Double.PositiveInfinity)
    state = state.copy(
      cycleBasis = state.cycleBasis + (simplex -> Chain(
        SortedMap(simplex -> 1.0)(using Ordering.by(filtrationValues.getOrElse(_, Double.PositiveInfinity)))
      )),
      cycleBirths = state.cycleBirths + (simplex -> birth)
    )
  } else {
    val dsimplex = Chain.boundary(simplex)(Ordering.by(filtrationValues.getOrElse(_, Double.PositiveInfinity)))
    val dimension = simplex.size - 1
    val (dsModB, reductionB) = Chain.reduceBy(dsimplex, state.boundaryBasis)
    if (dsModB.isZero) {
      // then ∂s is already a boundary, say the boundary of z
      // hence ∂(s-z) = ∂s - ∂z = 0
      // so s-z is a new cycle
      val cycle = Chain(
        SortedMap.from(reductionB.map((s, v) => (state.coboundaryLookup(s), v)))(using
          Ordering.by(filtrationValues.getOrElse(_, Double.PositiveInfinity))
        ) + (simplex -> 1.0)
      )
      val (cyclemodZ, _) = Chain.reduceBy(cycle, state.cycleBasis)
      val cycleBirth = filtrationValues.getOrElse(simplex, Double.PositiveInfinity)

      state = cyclemodZ.entries.headOption match {
        case Some((spx, _)) =>
          state.copy(cycleBasis = state.cycleBasis + (spx -> cyclemodZ), cycleBirths = state.cycleBirths + (spx -> cycleBirth))
        case None: Option[(Set[Int], Double)] => state
      }
    } else {
      // since ∂s is not a boundary already, it has to be a cycle
      // reducing mod cycles will tell us which cycle just became a boundary with the addition of simplex
      val (dsModZ, reductionZ) = Chain.reduceBy(dsModB, state.cycleBasis)
      assert(dsModZ.isZero)
      val dyingChainHead = reductionZ.keys
        .maxBy(filtrationValues.getOrElse(_, Double.PositiveInfinity))
      barcodes.append((dyingChainHead.size - 1, state.cycleBirths(dyingChainHead), filtrationValues(simplex)))
      state = state.copy(
        cycleBasis = state.cycleBasis - dyingChainHead,
        boundaryBasis = state.boundaryBasis + (dyingChainHead -> state.cycleBasis(dyingChainHead)),
        cycleBirths = state.cycleBirths - dyingChainHead,
        coboundaryLookup = state.coboundaryLookup + (dyingChainHead -> simplex)
      )
    }
  }
  def barcode: Seq[(Int, Double, Double)] = {
    simplexstream.foreach(step)
    state.cycleBirths
      .map { case (simplex, birth) =>
        ((simplex.size - 1, birth, Double.PositiveInfinity))
      }
      .toSeq
      .concat(barcodes)
  }
end NaivePersistentHomology
