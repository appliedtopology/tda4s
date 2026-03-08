package org.appliedtopology.tda4s

import Ordering.Implicits.*
import scala.annotation.tailrec
import scala.collection.immutable.SortedMap

case class Chain(entries: SortedMap[Set[Int], Double]):
  def scaleAdd(scale: Double, other: Chain): Chain =
    Chain(other.entries.foldLeft(entries) { case (result, (s, v)) =>
      entries.getOrElse(s, 0.0) match {
        case v0: Double =>
          v0 + scale * v match {
            case v1 if math.abs(v1) < 1e-10 => result.removed(s)
            case v1                         => result.updated(s, v1)
          }
      }
    })
  def isZero: Boolean = entries.values.forall(math.abs(_) < 1e-10)

object Chain:
  def boundary(simplex: Set[Int])(ordering: Ordering[Set[Int]]): Chain =
    Chain(
      SortedMap.from(
        simplex.toList.sorted
          .map(i => simplex - i)
          .zip(Iterator.iterate(1.0)(-_))
      )(using ordering)
    )
  @tailrec
  def reduceBy(chain: Chain, basis: Map[Set[Int], Chain], log: Map[Set[Int], Double] = Map.empty): (Chain, Map[Set[Int], Double]) =
    (for {
      (leadingSimplex, leadingValue) <- chain.entries.headOption
      basisChain <- basis.get(leadingSimplex)
      (basisSimplex, basisValue) <- basisChain.entries.headOption
    } yield (basisSimplex, -leadingValue / basisValue, chain.scaleAdd(-leadingValue / basisValue, basisChain))) match {
      case Some((basisSimplex, coeff, reducedChain)): Option[(Set[Int], Double, Chain)] =>
        reduceBy(reducedChain, basis, log.updated(basisSimplex, coeff))
      case None: Option[Chain] => (chain, log)
    }
