-- |
-- Module      : Data.VectorSpace.SIMD.FiniteSupportedSequence
-- Copyright   : (c) Justus Sagemüller 2017
-- License     : GPL v3
-- 
-- Maintainer  : (@) sagemueller $ geo.uni-koeln.de
-- Stability   : experimental
-- Portability : portable
-- 
{-# LANGUAGE TypeFamilies            #-}
{-# LANGUAGE TypeOperators           #-}
{-# LANGUAGE FlexibleInstances       #-}
{-# LANGUAGE FlexibleContexts        #-}
{-# LANGUAGE MultiParamTypeClasses   #-}
{-# LANGUAGE DefaultSignatures       #-}
{-# LANGUAGE ScopedTypeVariables     #-}
{-# LANGUAGE UnicodeSyntax           #-}
{-# LANGUAGE CPP                     #-}
{-# LANGUAGE PatternSynonyms         #-}
{-# LANGUAGE ViewPatterns            #-}

module Data.VectorSpace.SIMD.FiniteSupportedSequence where

import Math.DivisionAlgebra.SIMD

import Data.AdditiveGroup
import Data.VectorSpace
import Data.VectorSpace.SIMD

import qualified Data.Vector.Generic as Arr
import qualified Data.Vector.Generic.Mutable as MArr
import qualified Data.Vector as BArr
import qualified Data.Vector.Unboxed as UArr
import qualified Data.Vector.Fusion.Bundle as ABdl

import qualified Data.Primitive.SIMD as SIMD

import GHC.Exts (IsList(..))

trimTrailingZeroes :: (Num t, Eq t, Arr.Vector v t) => v t -> v t
trimTrailingZeroes v = Arr.force $ Arr.take nlnz v
 where nlnz = case ABdl.findIndex (/=0) $ Arr.streamR v of
         Just nrz -> Arr.length v - nrz
         Nothing  -> 0

class AdditiveGroup (FinSuppSeq t v) => PackSequence t v where
  data FinSuppSeq t v :: *
  toArray :: FinSuppSeq t v -> BArr.Vector v
  fromArray :: BArr.Vector v -> FinSuppSeq t v

instance PackSequence ℝ ℝ where
  data FinSuppSeq ℝ ℝ = ℝFinSuppSeqℝ Int (UArr.Vector ℝ)
  toArray (ℝFinSuppSeqℝ nLeadingZeroes v)
      = Arr.replicate nLeadingZeroes 0 Arr.++ Arr.convert v
  fromArray vb = ℝFinSuppSeqℝ (Arr.length leadingZeroes) $ trimTrailingZeroes vp
   where (leadingZeroes, vp) = UArr.span (==0) $ Arr.convert vb

instance AdditiveGroup (FinSuppSeq ℝ ℝ) where
  zeroV = ℝFinSuppSeqℝ 0 $ Arr.empty
  negateV (ℝFinSuppSeqℝ i₀ v) = ℝFinSuppSeqℝ i₀ $ Arr.map negate v
  ℝFinSuppSeqℝ i₀₀ v₀ ^+^ ℝFinSuppSeqℝ i₀₁ v₁ = uncurry ℝFinSuppSeqℝ
     $ addFinsuppSeqs (i₀₀, v₀) (i₀₁, v₁)

addFinsuppSeqs :: (Num n, UArr.Unbox n)
           => (Int, UArr.Vector n) -> (Int, UArr.Vector n) -> (Int, UArr.Vector n)
addFinsuppSeqs (i₀₀, v₀) (i₀₁, v₁) = (i₀s, vs)
   where vs = case (compare i₀₀ i₀₁, compare (i₀₀+Arr.length v₀) (i₀₁+Arr.length v₁)) of
               (EQ,EQ) -> Arr.zipWith (+) v₀ v₁
               (EQ,LT) -> Arr.zipWith (+) v₀ v₁ Arr.++ Arr.drop (Arr.length v₀) v₁
               (EQ,GT) -> Arr.zipWith (+) v₀ v₁ Arr.++ Arr.drop (Arr.length v₁) v₀
               (LT,LT) -> 
                   let δi₀ = i₀₁-i₀₀
                   in Arr.take δi₀ v₀
                    Arr.++ Arr.zipWith (+) (Arr.drop δi₀ v₀) v₁
                    Arr.++ Arr.drop (Arr.length v₀ - δi₀) v₁
               (GT,GT) -> 
                   let δi₀ = i₀₀-i₀₁
                   in Arr.take δi₀ v₁
                    Arr.++ Arr.zipWith (+) v₀ (Arr.drop δi₀ v₁)
                    Arr.++ Arr.drop (Arr.length v₁ - δi₀) v₀
               (LT,_) -> Arr.create $ do
                    zs <- Arr.thaw v₀
                    let δi₀ = i₀₁-i₀₀
                    Arr.imapM_ (\i x -> MArr.unsafeModify zs (+x) (i+δi₀)) v₁
                    return zs
               (GT,_) -> Arr.create $ do
                    zs <- Arr.thaw v₁
                    let δi₀ = i₀₀-i₀₁
                    Arr.imapM_ (\i x -> MArr.unsafeModify zs (+x) (i+δi₀)) v₀
                    return zs
         i₀s = min i₀₀ i₀₁
  
instance PackSequence t v => IsList (FinSuppSeq t v) where
  type Item (FinSuppSeq t v) = v
  fromList = fromArray . Arr.fromList
  toList = Arr.toList . toArray

instance (PackSequence t v, Show v) => Show (FinSuppSeq t v) where
  show = show . toList




instance PackSequence SIMD.DoubleX4 ℝ where
  data FinSuppSeq SIMD.DoubleX4 ℝ = ℝ⁴FinSuppSeqℝ Int (UArr.Vector SIMD.DoubleX4)
  toArray (ℝ⁴FinSuppSeqℝ nLeading4Zeroes v)
      = Arr.replicate (4*nLeading4Zeroes) 0
          Arr.++ (Arr.convert v >>= lrep . SIMD.unpackVector)
   where lrep (x,y,z,w) = Arr.fromList [x,y,z,w]
  fromArray vb = ℝ⁴FinSuppSeqℝ (Arr.length leading4Zeroes) $ trimTrailingZeroes vp
   where (leading4Zeroes, vp) = UArr.span (==0) vpacked
         vpacked = Arr.generate nby4
                   $ \i' -> let i = i'*4
                                at j | j<Arr.length vb  = Arr.unsafeIndex vb j
                                     | otherwise        = 0
                            in SIMD.packVector (at i, at (i+1), at (i+2), at (i+3))
         nby4 = (Arr.length vb + 3)`div`4

instance AdditiveGroup (FinSuppSeq SIMD.DoubleX4 ℝ) where
  zeroV = ℝ⁴FinSuppSeqℝ 0 $ Arr.empty
  negateV (ℝ⁴FinSuppSeqℝ i₀ v) = ℝ⁴FinSuppSeqℝ i₀ $ Arr.map negate v
  ℝ⁴FinSuppSeqℝ i₀₀ v₀ ^+^ ℝ⁴FinSuppSeqℝ i₀₁ v₁ = uncurry ℝ⁴FinSuppSeqℝ
     $ addFinsuppSeqs (i₀₀, v₀) (i₀₁, v₁)
