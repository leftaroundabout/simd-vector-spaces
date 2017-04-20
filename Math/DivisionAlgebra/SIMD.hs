-- |
-- Module      : Math.DivisionAlgebra.SIMD
-- Copyright   : (c) Justus Sagemüller 2017
-- License     : GPL v3
-- 
-- Maintainer  : (@) sagemueller $ geo.uni-koeln.de
-- Stability   : experimental
-- Portability : portable
-- 
{-# LANGUAGE TypeFamilies            #-}
{-# LANGUAGE FlexibleInstances       #-}
{-# LANGUAGE FlexibleContexts        #-}
{-# LANGUAGE ScopedTypeVariables     #-}
{-# LANGUAGE UnicodeSyntax           #-}

module Math.DivisionAlgebra.SIMD where

import Data.AffineSpace
import Data.VectorSpace
import Data.Basis

import qualified Data.Foldable as Foldable

import qualified Data.Primitive.SIMD as SIMD


data family Complex a

newtype instance Complex Double = ComplexDouble { getComplexFloat :: SIMD.DoubleX2 }

instance Num (Complex Double) where
  {-# INLINE fromInteger #-}
  fromInteger n = ComplexDouble $ SIMD.packVector (fromInteger n, 0)
  {-# INLINE (+) #-}
  ComplexDouble a + ComplexDouble b = ComplexDouble $ a + b
  {-# INLINE (-) #-}
  ComplexDouble a - ComplexDouble b = ComplexDouble $ a - b
  {-# INLINE (*) #-}
  ComplexDouble a * ComplexDouble b = ComplexDouble
      $ a * SIMD.broadcastVector rb + SIMD.packVector (-ia, ra) * SIMD.broadcastVector ib
   where (ra, ia) = SIMD.unpackVector a
         (rb, ib) = SIMD.unpackVector b
  {-# INLINE negate #-}
  negate (ComplexDouble a) = ComplexDouble $ negate a
  {-# INLINE abs #-}
  abs (ComplexDouble a) = ComplexDouble $ SIMD.packVector (sqrt . SIMD.sumVector $ a*a, 0)
  {-# INLINE signum #-}
  signum (ComplexDouble a) = ComplexDouble
                $ a / SIMD.broadcastVector (sqrt . SIMD.sumVector $ a*a)

instance AdditiveGroup (Complex Double) where
  {-# INLINE zeroV #-}
  zeroV = ComplexDouble SIMD.nullVector
  {-# INLINE (^+^) #-}
  (^+^) = (+)
  {-# INLINE negateV #-}
  negateV = negate

-- ^ Note that the scalars are complex, too!
instance VectorSpace (Complex Double) where
  type Scalar (Complex Double) = Complex Double
  {-# INLINE (*^) #-}
  (*^) = (*)
  
instance Show (Complex Double) where
  showsPrec p (ComplexDouble a) = showParen (p>6) $ shows r . ("+:"++) . shows i
   where (r,i) = SIMD.unpackVector a


infixl 6 +:
class VectorSpace c => DivisionAlgebra c where
  type RealPart c :: *
  type ImagPart c :: *
  (+:) :: RealPart c -> ImagPart c -> c
  realPart :: c -> RealPart c
  imagPart :: c -> ImagPart c

instance DivisionAlgebra Double where
  type RealPart Double = Double
  type ImagPart Double = ()
  x+:() = x
  realPart = id
  imagPart = const ()
instance DivisionAlgebra (Complex Double) where
  type RealPart (Complex Double) = Double
  type ImagPart (Complex Double) = Double
  x+:y = ComplexDouble $ SIMD.packVector (x,y)
  realPart (ComplexDouble a) = r
   where (r, _) = SIMD.unpackVector a
  imagPart (ComplexDouble a) = i
   where (_, i) = SIMD.unpackVector a
