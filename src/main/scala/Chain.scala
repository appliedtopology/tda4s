package org.appliedtopology.tda4s

import Ordering.Implicits.*
import scala.annotation.tailrec
import scala.collection.immutable.SortedMap
import scala.collection.{Map, Set}

/** Chains carry information about how simplexes are connected. They are formal linear combinations of simplexes with coefficients drawn
  * over a field. Our implementation is parameterized by a type T which is then typecasted to the Field type.
  *
  * @tparam T
  *   The underlying type for the field coefficients. Should be double for the reals and int for a finite field
  *
  * @param simplexCoefficients
  *   Chains are implemented as conversions of a SortedMap from the simplex to the corresponding coefficient. Sorting should occur by
  *   filtration value in the simplex stream.
  */
class Chain[T](using val field: Field[T])(val simplexCoefficients: SortedMap[Simplex, field.F]):
  import field.*

  /** Scales all the coefficients in the chain by a constant factor
    *
    * @param scalingFactor
    *   the field element which is used to multiply each coefficient
    */
  def scale(scalingFactor: F): Chain[T] = Chain(simplexCoefficients.collect[Simplex, F] { case (k, v) =>
    (k, v * scalingFactor)
  }(using simplexCoefficients.ordering))

  /** Constructs a new chain by scaling the original (self) chain and adding it to the other chain
    *
    * @param scalingFactor
    *   the field element which is used to multiply each coefficient
    * @param other
    *   the other chain used for addition coefficient-wise
    */
  def scaleAdd(scalingFactor: F, other: Chain[T]): Chain[T] =
    Chain(other.simplexCoefficients.foldLeft(simplexCoefficients) { case (result, (s, v)) =>
      result.getOrElse(s, zero) match {
        case v0: F =>
          v match {
            case v1: F =>
              v0 + scalingFactor * v1 match {
                case v2 if field.isZero(v2) => result.removed(s)
                case v2                     => result.updated(s, v2)
              }
          }
      }
    })
  def isZero: Boolean = simplexCoefficients.values.forall(field.isZero(_))
  def leading: Option[(Simplex, F)] = simplexCoefficients.lastOption

//Companion object for the Chain Class.
//We utilize it to create a class factory, particularly to construct chains from simplexes.
object Chain:
  def from[T: Field as field](simplex: Simplex)(ordering: Ordering[Simplex]): Chain[T] = Chain(
    SortedMap((simplex, field.one))(using ordering)
  )

  /** The from method allows for the construction of the Chain object from the components
    *
    * @tparam Field
    *   [T] type bound for our field, requires us to have an instantiation of F[T]
    * @param terms
    *   Collection of Simplex and field element to store in the chain map
    * @param ordering
    *   We need a way to define an ordering on simplexes for the persistent homology algorithm to work
    */
  def from[T: Field as field](terms: (Simplex, field.F)*)(using ordering: Ordering[Simplex]): Chain[T] =
    Chain(using field)(SortedMap.from(terms)(using ordering))
  def reduce[T: Field as outerField](chain: Chain[T], basis: Map[Simplex, Chain[T]]): (Chain[T], Map[Simplex, outerField.F]) =
    import outerField._
    @tailrec
    def reduceBy(
        innerChain: Chain[T],
        basis: Map[Simplex, Chain[T]],
        reductionHistory: Map[Simplex, F]
    )(using innerField: Field[T]): (Chain[T], Map[Simplex, F]) = innerChain.field match {
      case outerField =>
        innerChain.simplexCoefficients.find(k => basis.contains(k._1)) match {
          case None => (innerChain, reductionHistory)
          case Some((simplex, chainCoeff: F)) =>
            basis.get(simplex) match {
              case None => (innerChain, reductionHistory)
              case Some(pivotChain) =>
                val pivotCoeff: F = pivotChain.simplexCoefficients.getOrElse(simplex, zero).asInstanceOf[F]
                // We use a type conversion here since Scala was raising the issue that chains might draw coefficients from different fields
                val resultCoeff: F = chainCoeff / pivotCoeff
                reduceBy(
                  innerChain.scaleAdd(-resultCoeff.asInstanceOf[innerChain.field.F], pivotChain),
                  basis,
                  reductionHistory + (simplex -> resultCoeff)
                )(using outerField)
            }
        }
    }
    reduceBy(chain, basis, Map.empty)(using outerField)
