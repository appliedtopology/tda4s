package org.appliedtopology.tda4s

import language.experimental.modularity
import org.scalacheck.Properties
import org.scalacheck.Prop.{forAll, propBoolean}
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.propspec.AnyPropSpec
import org.scalatestplus.scalacheck.ScalaCheckPropertyChecks
import org.scalatest.matchers.should.Matchers
import org.scalatest.matchers.should.Matchers.shouldEqual

class PersistentHomologySpec extends AnyFlatSpec:
  "The tetrahedron" should "have barcode" in {
    val simplexstream = Seq(
      Set(0),
      Set(1),
      Set(2),
      Set(3),
      Set(0, 1),
      Set(0, 2),
      Set(0, 3),
      Set(1, 2),
      Set(1, 3),
      Set(2, 3),
      Set(0, 1, 2),
      Set(0, 1, 3),
      Set(0, 2, 3),
      Set(1, 2, 3)
    )
    val filtrationValues = Map(
      Set(0) -> 0.0,
      Set(1) -> 1.0,
      Set(2) -> 2.0,
      Set(3) -> 3.0,
      Set(0, 1) -> 4.0,
      Set(0, 2) -> 5.0,
      Set(0, 3) -> 6.0,
      Set(1, 2) -> 7.0,
      Set(1, 3) -> 8.0,
      Set(2, 3) -> 9.0,
      Set(0, 1, 2) -> 10.0,
      Set(0, 1, 3) -> 11.0,
      Set(0, 2, 3) -> 12.0,
      Set(1, 2, 3) -> 13.0
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
    PersistentHomology.persistentHomology(simplexstream, filtrationValues).sorted shouldEqual expectedBarcode.sorted
  }
