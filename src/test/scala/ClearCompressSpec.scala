package org.appliedtopology.tda4s

import org.scalacheck.Properties
import org.scalacheck.Prop.{forAll, propBoolean}
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.propspec.AnyPropSpec
import org.scalatestplus.scalacheck.{Checkers, ScalaCheckPropertyChecks}
import org.scalatest.matchers.should.Matchers
import org.scalatest.matchers.should.Matchers.shouldEqual

import org.apache.commons.rng.simple.RandomSource

import scala.math.{cos, sin}
import scala.collection.Set

class ClearCompressSpec extends AnyFlatSpec:
  "The tetrahedron" should "have barcode" in {
    val simplexstream = Seq(
      Simplex(0),
      Simplex(1),
      Simplex(2),
      Simplex(3),
      Simplex(0, 1),
      Simplex(0, 2),
      Simplex(0, 3),
      Simplex(1, 2),
      Simplex(1, 3),
      Simplex(2, 3),
      Simplex(0, 1, 2),
      Simplex(0, 1, 3),
      Simplex(0, 2, 3),
      Simplex(1, 2, 3)
    ).iterator
    val filtrationValues: PartialFunction[Simplex, Double] = Map(
      Simplex(0) -> 0.0,
      Simplex(1) -> 1.0,
      Simplex(2) -> 2.0,
      Simplex(3) -> 3.0,
      Simplex(0, 1) -> 4.0,
      Simplex(0, 2) -> 5.0,
      Simplex(0, 3) -> 6.0,
      Simplex(1, 2) -> 7.0,
      Simplex(1, 3) -> 8.0,
      Simplex(2, 3) -> 9.0,
      Simplex(0, 1, 2) -> 10.0,
      Simplex(0, 1, 3) -> 11.0,
      Simplex(0, 2, 3) -> 12.0,
      Simplex(1, 2, 3) -> 13.0
    )
    val expectedBarcode = Seq(
      (0, 0.0, Double.PositiveInfinity),
      (0, 1.0, 4.0),
      (0, 2.0, 5.0),
      (0, 3.0, 6.0),
      (1, 7.0, 10.0),
      (1, 8.0, 11.0),
      (1, 9.0, 12.0),
      (2, 13.0, Double.PositiveInfinity)
    )
    given DoubleIsField(1e-15)
    ClearCompress.persistentHomology(simplexstream, filtrationValues).sorted shouldEqual expectedBarcode.sorted
  }

val rngCC = RandomSource.MT.create()

class ClearCompressVietorisRipsValidation extends AnyFlatSpec with Checkers with Matchers:
  "the circle" should "have barcode" in {
    given FiniteFieldIsField(17)
    val D = 3
    val N = 10
    val aa = (0 until N).map(_ / N.toDouble)
    val xy = aa.map(a => Array(cos(2 * 3.14 * a) + rngCC.nextDouble(-0.05, 0.05), sin(2 * 3.14 * a) + rngCC.nextDouble(-0.05, 0.05))).toSeq
    val metricSpace = VectorMetricSpace("euclidean", xy)
    val vrStream = NaiveVietorisRips(metricSpace)
    val barcode = ClearCompress.persistentHomology(vrStream.simplices().filter(_.dimension <= D), vrStream.filtrationValues)
    val salientBars = barcode.filter((bar: (Int, Double, Double)) => bar._1 < D)
    val numBars = salientBars
      .filter((bar: (Int, Double, Double)) => (bar._3 - bar._2) > 0.75)
      .size
    numBars shouldEqual 2
  }

class ClearCompressAlphashapeValidation extends AnyFlatSpec with Checkers with Matchers:
  "the circle" should "have barcode" in {
    given FiniteFieldIsField(17)
    val D = 2
    val N = 10
    val aa = (0 until N).map(_ / N.toDouble)
    val xy = aa.map(a => Array(cos(2 * 3.14 * a) + rngCC.nextDouble(-0.05, 0.05), sin(2 * 3.14 * a) + rngCC.nextDouble(-0.05, 0.05))).toSeq
    val alpha = Alpha(xy)
    val barcode = ClearCompress.persistentHomology(alpha.simplices().filter(_.dimension <= D), alpha.filtrationValues)
    val salientBars = barcode.filter((bar: (Int, Double, Double)) => bar._1 < D)
    val numBars = salientBars
      .filter((bar: (Int, Double, Double)) => (bar._3 - bar._2) > 0.5)
      .size
    numBars shouldEqual 2
  }

class ClearCompressKleinBottleValidation extends AnyFlatSpec with Checkers with Matchers:
  "the Klein bottle" should "have different barcodes" in {
    val topSimplices: Seq[Simplex] = Seq(
      Simplex(2, 6, 7),
      Simplex(0, 6, 7),
      Simplex(2, 5, 7),
      Simplex(0, 5, 7),
      Simplex(4, 5, 6),
      Simplex(3, 5, 6),
      Simplex(0, 4, 6),
      Simplex(2, 3, 6),
      Simplex(1, 4, 5),
      Simplex(0, 3, 5),
      Simplex(1, 2, 5),
      Simplex(2, 3, 4),
      Simplex(1, 3, 4),
      Simplex(0, 2, 4),
      Simplex(0, 1, 3),
      Simplex(0, 1, 2)
    )
    val simplices = topSimplices.toSet
      .flatMap((s: Simplex) => s.vertices.subsets())
      .map(Simplex(_))
      .filter(_.dimension >= 0)
      .toSeq
      .sorted(using Ordering.by((s: Simplex) => s.vertices)(using Ordering.Implicits.sortedSetOrdering))
      .sortBy(_.dimension)
    println(simplices)
    val filtrationValues: PartialFunction[Simplex, Double] = { case _ => 0.0 }
    val barcode0 = {
      given DoubleIsField(1e-15)
      ClearCompress.persistentHomology(simplices.iterator, filtrationValues)
    }
    val barcode2 = {
      given FiniteFieldIsField(2)
      ClearCompress.persistentHomology(simplices.iterator, filtrationValues)
    }
    val barcode3 = {
      given FiniteFieldIsField(3)
      ClearCompress.persistentHomology(simplices.iterator, filtrationValues)
    }
    println(barcode0)
    println(barcode2)
    println(barcode3)
    assert(barcode0 == barcode3)
    assert(barcode0 != barcode2)
    assert(Set(barcode0*) == Set((0, 0.0, Double.PositiveInfinity), (1, 0.0, Double.PositiveInfinity)))
    assert(
      Set(barcode2*) == Set(
        (0, 0.0, Double.PositiveInfinity),
        (1, 0.0, Double.PositiveInfinity),
        (1, 0.0, Double.PositiveInfinity),
        (2, 0.0, Double.PositiveInfinity)
      )
    )
    assert(Set(barcode3*) == Set((0, 0.0, Double.PositiveInfinity), (1, 0.0, Double.PositiveInfinity)))
  }
