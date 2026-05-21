package org.appliedtopology.tda4s

import collection.mutable

// ChunkedClearCompressWorker
//
// Implements the chunked clear and compress optimization for persistent homology
// boundary matrix reduction. The filtration is divided into chunks, and reduction
// proceeds in three phases.
//
// First, localReduction() performs standard matrix reduction
// within each chunk, using clearing to zero out positive columns
//
// Second, compress() eliminates entries from unpaired
// columns that cannot affect the remaining global persistence pairs.
//
// Finally, reduceGlobal() performs reduction using again the clear optimization on the remaining global matrix, which
// is now substantially smaller due to the earlier local and compression phases.

object ClearCompress:
  def persistentHomology[T: Field](
      simplexStream: Iterator[Simplex],
      filtrationValues: PartialFunction[Simplex, Double]
  ): Seq[(Int, Double, Double)] =
    val simplicesSeq = simplexStream.toIndexedSeq
    val chunkSize = math.sqrt(simplicesSeq.size).toInt.max(1)
    val worker = ChunkedClearCompressWorker[T](simplicesSeq.iterator, filtrationValues, chunkSize)
    val boundaries = mutable.Map.empty[Simplex, Chain[T]]
    worker.localReduction(boundaries)
    worker.markActiveEntries()
    worker.compress(boundaries)
    worker.reduceGlobal(boundaries)
    val getVal = filtrationValues.applyOrElse(_, _ => Double.PositiveInfinity)
    val finite = worker.pairOf.toSeq
      .filter((birth, _) => !worker.negatives.contains(birth))
      .map((birth, death) => (birth.dimension, getVal(birth), getVal(death)))
      .filter { case (_, b, d) => math.abs(d - b) > 1e-15 }
    val infinite = worker.allCells
      .filter(s => !worker.paired.contains(s))
      .map(s => (s.dimension, getVal(s), Double.PositiveInfinity))
    (finite ++ infinite).toSeq
end ClearCompress

class ChunkedClearCompressWorker[CoefficientT: Field as fr](
    simplexIter: Iterator[Simplex],
    val filtrationValues: PartialFunction[Simplex, Double],
    chunkSize: Int
):
  import fr.*

  def filtrationValue(simplex: Simplex): Double =
    filtrationValues.applyOrElse(simplex, _ => Double.PositiveInfinity)

  val filtrationOrdering: Ordering[Simplex] =
    Ordering.by(filtrationValue).orElseBy(_.vertices)(Ordering.Implicits.sortedSetOrdering)

  given Ordering[Simplex] = filtrationOrdering

  def emptyChain: Chain[CoefficientT] = Chain.from[CoefficientT]()

  def boundaryOf(sigma: Simplex): Chain[CoefficientT] =
    if sigma.dimension <= 0 then emptyChain
    else sigma.boundary[CoefficientT](filtrationOrdering)

  private def columnReduce(col: Chain[CoefficientT], ratio: F, pivot: Chain[CoefficientT]): Chain[CoefficientT] =
    col.scaleAdd(neg(ratio).asInstanceOf[col.field.F], pivot)

  private def leadCoeff(c: Chain[CoefficientT]): F =
    c.leading.get._2.asInstanceOf[F]

  val paired: mutable.Set[Simplex] = mutable.Set.empty
  val cleared: mutable.Set[Simplex] = mutable.Set.empty
  val negatives: mutable.Set[Simplex] = mutable.Set.empty
  val pairOf: mutable.Map[Simplex, Simplex] = mutable.Map.empty

  // Constructing chunks
  val allCells: IndexedSeq[Simplex] = simplexIter.toIndexedSeq.sorted(using filtrationOrdering)

  val chunks: IndexedSeq[IndexedSeq[Simplex]] =
    allCells.grouped(chunkSize).toIndexedSeq

  val chunkIndex: Map[Simplex, Int] =
    chunks.zipWithIndex.flatMap { (chunk, idx) =>
      chunk.map(cell => cell -> idx)
    }.toMap

  // Tracking reduced state of each column across phases.
  val columnState: mutable.Map[Simplex, Chain[CoefficientT]] = mutable.Map.empty

  def isLocal(sigma: Simplex, i: Int): Boolean =
    chunkIndex.get(sigma).exists(j => math.abs(j - i) <= 1)

  def isGlobal(sigma: Simplex, i: Int): Boolean = !isLocal(sigma, i)

  // Algorithm 2: LOCAL REDUCTION

  def localReduction(
      boundaries: mutable.Map[Simplex, Chain[CoefficientT]]
  ): Unit =

    val maxDim = allCells.map(_.dimension).maxOption.getOrElse(0)
    val numChunks = chunks.length

    // compute position of cells for fast lookup
    val cellPos: Map[Simplex, Int] = allCells.zipWithIndex.toMap

    // returns the position of the last cell
    def windowBoundary(k: Int): Int =
      if k <= 0 then -1
      else chunks(k - 1).lastOption.flatMap(cellPos.get).getOrElse(-1)

    // process cells from highest to lowest dimension,
    for dim <- maxDim.to(0).by(-1) do
      // pass=1 handles cells within a single chunk
      // pass=2 handles cells whose pivot exist in previous chunk
      for pass <- 1 to 2 do
        for chunk <- pass to numChunks do
          val chunkCells = chunks(chunk - 1)
          val t_br = windowBoundary(chunk - pass)

          for sigma <- chunkCells if sigma.dimension == dim do
            if !paired.contains(sigma) then

              var Rj: Chain[CoefficientT] =
                columnState.getOrElse(sigma, boundaryOf(sigma))

              // standard column reduction
              var reducing = true
              while reducing do
                Rj.leading match
                  case None =>
                    reducing = false
                  case Some((pivot, _)) =>
                    val pivotPos = cellPos.getOrElse(pivot, -1)
                    if pivotPos > t_br && boundaries.contains(pivot) then
                      val pivotChain = boundaries(pivot)
                      val ratio = div(leadCoeff(Rj), leadCoeff(pivotChain))
                      Rj = columnReduce(Rj, ratio, pivotChain)
                    else reducing = false

              columnState(sigma) = Rj

              // recording persistent pairs found locally
              if !Rj.isZero then
                Rj.leading match
                  case Some((i, _)) =>
                    val iPos = cellPos.getOrElse(i, -1)
                    if iPos > t_br then
                      boundaries(i) = Rj
                      cleared.add(i)
                      columnState(i) = emptyChain
                      paired.add(i)
                      paired.add(sigma)
                      negatives.add(sigma)
                      pairOf(i) = sigma
                      pairOf(sigma) = i
                  case None => ()

  // Algorithm 3: MARK ACTIVE ENTRIES

  val active: mutable.Set[Simplex] = mutable.Set.empty
  val inactive: mutable.Set[Simplex] = mutable.Set.empty

  def markColumn(col: Simplex): Boolean =
    import scala.util.boundary
    import scala.util.boundary.break

    if active.contains(col) then return true
    if inactive.contains(col) then return false

    val Rk: Chain[CoefficientT] =
      columnState.getOrElse(col, boundaryOf(col))

    // scan entries.If any entry is unpaired, or paired with an active column,
    // then this column is active and must be kept.
    val foundActive: Boolean = boundary:
      for (ell, _) <- Rk.simplexCoefficients do
        if !paired.contains(ell) then break(true)
        else if !negatives.contains(ell) then
          pairOf.get(ell) match
            case Some(j) if j != col && markColumn(j) =>
              break(true)
            case _ => ()
      false

    if foundActive then
      active.add(col)
      true
    else
      inactive.add(col)
      false

  // run markColumn on every unpaired, uncleared cell
  def markActiveEntries(): Unit =
    for sigma <- allCells if !paired.contains(sigma) && !cleared.contains(sigma) do markColumn(sigma)

  // Algorithm 4: COMPRESS

  def compress(
      boundaries: mutable.Map[Simplex, Chain[CoefficientT]]
  ): Unit =

    for sigma <- allCells if !paired.contains(sigma) && !cleared.contains(sigma) do

      var Rj: Chain[CoefficientT] =
        columnState.getOrElse(sigma, boundaryOf(sigma))

      // remove entries that are not connected to any global columns.
      // inactive pivot columns can be dropped
      Rj = Rj.simplexCoefficients.foldLeft(Rj) { case (acc, (ell, coeff)) =>
        val c = coeff.asInstanceOf[F]
        if negatives.contains(ell) then acc
        else if boundaries.contains(ell) && active.contains(ell) then
          val pivotCol = boundaries(ell)
          val ratio = div(c, leadCoeff(pivotCol))
          columnReduce(acc, ratio, pivotCol)
        else acc
      }

      columnState(sigma) = Rj

  // Algorithm 5, REDUCE GLOBAL SUBMATRIX

  def reduceGlobal(
      boundaries: mutable.Map[Simplex, Chain[CoefficientT]]
  ): Unit =

    val maxDim = allCells.map(_.dimension).maxOption.getOrElse(0)
    // only process cells that weren't already paired
    for dim <- maxDim.to(0).by(-1) do
      for
        sigma <- allCells
        if !paired.contains(sigma) && !cleared.contains(sigma) && sigma.dimension == dim
      do
        var Rj: Chain[CoefficientT] =
          columnState.getOrElse(sigma, boundaryOf(sigma))

        // same standard reduction as local phase but now over the full global matrix
        var reducing = true
        while reducing do
          Rj.leading match
            case None =>
              reducing = false
            case Some((pivot, _)) =>
              if boundaries.contains(pivot) then
                val pivotChain = boundaries(pivot)
                val ratio = div(leadCoeff(Rj), leadCoeff(pivotChain))
                Rj = columnReduce(Rj, ratio, pivotChain)
              else reducing = false

        // record persistent pairs found in remaining global matrix
        if !Rj.isZero then
          Rj.leading match
            case Some((i, _)) =>
              boundaries(i) = Rj
              cleared.add(i)
              columnState(i) = emptyChain
              paired.add(i)
              paired.add(sigma)
              negatives.add(sigma)
              pairOf(i) = sigma
              pairOf(sigma) = i
            case None => ()
