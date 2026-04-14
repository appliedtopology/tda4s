package org.appliedtopology.tda4s

import scala.collection.SortedSet

/** Case Class for Simplex Object. Equipped with set difference (-), union (+), dimension, and a boundary map
  *
  * @param vertices
  *   Simplexes are constructed from a sorted set of ints
  */
case class Simplex(vertices: SortedSet[Int]):
  def +(vertex: Int): Simplex = Simplex(vertices.union(Set(vertex)))
  def -(vertex: Int): Simplex = Simplex(vertices.diff(SortedSet(vertex)))
  def dimension: Int = vertices.size - 1

  /** The boundary of a k simplex is a chain of (k-1) simplexes, and these can be used to form a chain complex
    *
    * @tparam Field
    *   [T] type bound for our field, requires us to have an instantiation of F[T]
    * @param ordering
    *   We need a way to define an ordering on simplexes for the sorted maps that store chains
    */
  def boundary[T: Field as field](ordering: Ordering[Simplex]): Chain[T] = {
    val faces = vertices.toSeq.map(i => this - i)
    val coefficients: Iterator[field.F] = Iterator.iterate(field.one)(field.neg)
    val faceCoefficientPairs: Seq[(Simplex, field.F)] = faces.zip(coefficients).toSeq
    Chain.from(using field)(faceCoefficientPairs*)(using ordering)
  }

/** Companion object for the Simplex Class, used for several static methods relating to Simplexes
 * Particularly it's leveraged as a class factory to construct Simplexes from iterable collections
 * and destruct them as necessary
 */
object Simplex:
  def apply(vertices: Int*): Simplex = Simplex(SortedSet(vertices*))
  def from[T <: IterableOnce[Int]](vertices: T): Simplex = Simplex(SortedSet(vertices.iterator.toSeq*))
  def unapply(simplex: Simplex): Option[SortedSet[Int]] = Some(simplex.vertices)
  def unapplySeq(simplex: Simplex): Option[Seq[Int]] = Some(simplex.vertices.toSeq)
