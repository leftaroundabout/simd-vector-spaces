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

import Data.Int

import GHC.Exts (IsList(..))

import Test.QuickCheck.Arbitrary

trimTrailingZeroes :: (Num t, Eq t, Arr.Vector v t) => v t -> v t
trimTrailingZeroes v = Arr.force $ Arr.take nlnz v
 where nlnz = case ABdl.findIndex (/=0) $ Arr.streamR v of
         Just nrz -> Arr.length v - nrz
         Nothing  -> 0

class AdditiveGroup (FinSuppSeq t v) => PackSequence t v where
  data FinSuppSeq t v :: *
  toArray :: FinSuppSeq t v -> BArr.Vector v
  fromArray :: BArr.Vector v -> FinSuppSeq t v
  reoptimiseSequence :: FinSuppSeq t v -> FinSuppSeq t v
  reoptimiseSequence = id
  localBitDepth :: Functor p => p (t,v) -> Int

instance PackSequence ℝ ℝ where
  data FinSuppSeq ℝ ℝ = ℝFinSuppSeqℝ Int (UArr.Vector ℝ)
  toArray (ℝFinSuppSeqℝ nLeadingZeroes v)
      = Arr.replicate nLeadingZeroes 0 Arr.++ Arr.convert v
  fromArray vb = ℝFinSuppSeqℝ (Arr.length leadingZeroes) $ trimTrailingZeroes vp
   where (leadingZeroes, vp) = UArr.span (==0) $ Arr.convert vb
  localBitDepth _ = 52

instance AdditiveGroup (FinSuppSeq ℝ ℝ) where
  zeroV = ℝFinSuppSeqℝ 0 $ Arr.empty
  negateV (ℝFinSuppSeqℝ i₀ v) = ℝFinSuppSeqℝ i₀ $ Arr.map negate v
  ℝFinSuppSeqℝ i₀₀ v₀ ^+^ ℝFinSuppSeqℝ i₀₁ v₁ = uncurry ℝFinSuppSeqℝ
     $ addFinsuppSeqs (i₀₀, v₀) (i₀₁, v₁)

instance VectorSpace (FinSuppSeq ℝ ℝ) where
  type Scalar (FinSuppSeq ℝ ℝ) = ℝ
  μ *^ ℝFinSuppSeqℝ i₀ v
      = ℝFinSuppSeqℝ i₀ $ Arr.map (* μ) v

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
  localBitDepth _ = 52

instance AdditiveGroup (FinSuppSeq SIMD.DoubleX4 ℝ) where
  zeroV = ℝ⁴FinSuppSeqℝ 0 $ Arr.empty
  negateV (ℝ⁴FinSuppSeqℝ i₀ v) = ℝ⁴FinSuppSeqℝ i₀ $ Arr.map negate v
  ℝ⁴FinSuppSeqℝ i₀₀ v₀ ^+^ ℝ⁴FinSuppSeqℝ i₀₁ v₁ = uncurry ℝ⁴FinSuppSeqℝ
     $ addFinsuppSeqs (i₀₀, v₀) (i₀₁, v₁)

instance VectorSpace (FinSuppSeq SIMD.DoubleX4 ℝ) where
  type Scalar (FinSuppSeq SIMD.DoubleX4 ℝ) = ℝ
  μ *^ (ℝ⁴FinSuppSeqℝ i₀ v)
      = ℝ⁴FinSuppSeqℝ i₀ $ Arr.map (* SIMD.broadcastVector μ) v


type Exponent = Int

instance PackSequence SIMD.Int8X16 ℝ where
  data FinSuppSeq SIMD.Int8X16 ℝ
      = ℤ⁸'¹⁶FinSuppSeqℝ Int Exponent Int8 (UArr.Vector SIMD.Int8X16)
  toArray (ℤ⁸'¹⁶FinSuppSeqℝ nLeading16Zeroes magn _ v)
      = Arr.replicate (16*nLeading16Zeroes) 0
          Arr.++ (Arr.convert v >>= lrep . SIMD.unpackVector)
   where lrep (α,β,γ,δ,ε,ζ,η,θ,ι,κ,λ,μ,ν,ξ,ο,π)
             = Arr.map ((*2^^magn) . fromIntegral)
                 $ Arr.fromList [α,β,γ,δ,ε,ζ,η,θ,ι,κ,λ,μ,ν,ξ,ο,π]
  fromArray vb = ℤ⁸'¹⁶FinSuppSeqℝ (Arr.length leading16Zeroes) magn absMax
                   $ trimTrailingZeroes vp
   where (leading16Zeroes, vp) = UArr.span (==0) vpacked
         vpacked = Arr.generate nby4
                   $ \i' -> let i = i'*16
                                at j | j<Arr.length vb 
                                         = round $ Arr.unsafeIndex vb j / 2^^magn
                                     | otherwise        = 0
                            in SIMD.packVector
                                ( at i     , at (i+1),  at (i+2),  at (i+3)
                                , at (i+4) , at (i+5),  at (i+6),  at (i+7)
                                , at (i+8) , at (i+9),  at (i+10), at (i+11)
                                , at (i+12), at (i+13), at (i+14), at (i+15) )
         absMax = round $ origMax / 2^^magn
         origMax = Arr.maximum $ Arr.map abs vb
         bitDepth = fromIntegral $ localBitDepth [(Arr.head vpacked, origMax)]
         magn = floor $ logBase 2 (max tinyVal origMax) - bitDepth
         nby4 = (Arr.length vb + 3)`div`4
  reoptimiseSequence (ℤ⁸'¹⁶FinSuppSeqℝ nlz magn _ v)
                    = ℤ⁸'¹⁶FinSuppSeqℝ (nlz+nExtralz) magn' absMax v'
   where absMax = UArr.foldl' (\m -> max m . SIMD.foldVector max) 0 $ UArr.map abs v
         μ = max 0 $ round (logBase 2 (max tinyVal $ fromIntegral absMax) - bitDepth)
         bitDepth = fromIntegral $ localBitDepth [(Arr.head v, 0::ℝ)]
         nExtralz = maybe 0 id $ UArr.findIndex (/=0) v'₀
         magn' = undefined
         v'₀ | μ>0        = UArr.map (`SIMD.quotVector` SIMD.broadcastVector (2^μ)) v
             | otherwise  = v
         v' = trimTrailingZeroes $ UArr.drop nExtralz v'₀
  localBitDepth _ = 4

instance AdditiveGroup (FinSuppSeq SIMD.Int8X16 ℝ) where
  zeroV = ℤ⁸'¹⁶FinSuppSeqℝ 0 0 0 $ Arr.empty
  negateV (ℤ⁸'¹⁶FinSuppSeqℝ i₀ magn absMax v)
      = ℤ⁸'¹⁶FinSuppSeqℝ i₀ magn absMax $ Arr.map negate v
  ℤ⁸'¹⁶FinSuppSeqℝ i₀₀ magn₀ absMax₀ v₀ ^+^ ℤ⁸'¹⁶FinSuppSeqℝ i₀₁ magn₁ absMax₁ v₁
      = (if needsRenormalisation then reoptimiseSequence else id)
          $ ℤ⁸'¹⁶FinSuppSeqℝ i₀s magns absMaxs vs
   where (i₀s, vs) = addFinsuppSeqs (i₀₀, v₀') (i₀₁, v₁')
         [v₀', v₁'] = [ (if μ>0
                          then UArr.map (`SIMD.quotVector` SIMD.broadcastVector (2^μ))
                          else id) v
                      | (v,μ) <- [(v₀,μ₀),(v₁,μ₁)] ]
         absMaxs :: Int8
         absMaxs = sum [ absMax `div` (2^μ)
                       | (absMax,μ) <- [(absMax₀,μ₀), (absMax₁,μ₁)] ]
         μ₀, μ₁ :: Int
         ((μ₀, μ₁), (magns, needsRenormalisation))
           = case compare magn₀ magn₁ of
              LT | absMax₀ < maxBound`div`2
                    -> ((magn₁ - magn₀, 0), (magn₁, False))
              GT | absMax₁ < maxBound`div`2
                    -> ((0, magn₀ - magn₁), (magn₀, False))
              EQ | absMax₀ < maxBound`div`2 && absMax₁ < maxBound`div`2
                    -> ((0, 0), (magn₁, False))
              _  | maxMagn <- max magn₀ magn₁
                    -> ( (1+maxMagn-magn₀, 1+maxMagn-magn₁)
                       , (maxMagn-1, True) )




tinyVal :: Double
tinyVal = 3e-324 



instance Arbitrary (FinSuppSeq ℝ ℝ) where arbitrary = fmap fromList arbitrary
instance Arbitrary (FinSuppSeq SIMD.DoubleX4 ℝ) where arbitrary = fmap fromList arbitrary
instance Arbitrary (FinSuppSeq SIMD.Int8X16 ℝ) where arbitrary = fmap fromList arbitrary
