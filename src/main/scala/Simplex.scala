package org.appliedtopology.tda4s

import scala.collection.SortedSet

case class Simplex(vertices: SortedSet[Int]):
  def +(vertex: Int): Simplex = Simplex(vertices + vertex)
  def -(vertex: Int): Simplex = Simplex(vertices.diff(SortedSet(vertex)))
  def dimension: Int = vertices.size - 1
  def boundary[T: Field as field](ordering: Ordering[Simplex]): Chain[T] = {
    val faces = vertices.toSeq.map(i => this - i)
    val coefficients: Iterator[field.F] = Iterator.iterate(field.one)(field.neg)
    val faceCoefficientPairs: Seq[(Simplex, field.F)] = faces.zip(coefficients).toSeq
    Chain.from(using field)(faceCoefficientPairs*)(using ordering)
  }

object Simplex:
  def apply(vertices: Int*): Simplex = Simplex(SortedSet(vertices*))
  def from[T <: IterableOnce[Int]](vertices: T): Simplex = Simplex(SortedSet(vertices.iterator.toSeq*))
  def unapply(simplex: Simplex): Option[SortedSet[Int]] = Some(simplex.vertices)
  def unapplySeq(simplex: Simplex): Option[Seq[Int]] = Some(simplex.vertices.toSeq)
