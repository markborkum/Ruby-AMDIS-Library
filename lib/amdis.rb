###
# This module exports data types and parsers for Automatic Mass Spectral
# Deconvolution and Identification System (AMDIS) documents.
#
# Here is an example:
#
#   require 'amdis'
#   require 'open-uri'
#   
#   doc = AMDIS::MSL(open('http://www.example.com/path/to/msl/document'))
#   
#   doc.records.each do |record|
#     puts record
#   end
#
module AMDIS
  ###
  # Parse MSL.  Convenience method for `AMDIS::MSL::Document.parse`.
  #
  # @param object [Object]
  # @yieldparam doc [AMDIS::MSL::Document]
  # @return [AMDIS::MSL::Document]
  def self.MSL(object, &block)
    AMDIS::MSL::Document.parse(object, &block)
  end

  ###
  # This module exports data types and parsers for MSL documents.
  module MSL
    ###
    # Parse MSL.  Convenience method for `AMDIS::MSL::Document.parse`.
    #
    # @param object [Object]
    # @yieldparam doc [AMDIS::MSL::Document]
    # @return [AMDIS::MSL::Document]
    def self.parse(object, &block)
      AMDIS::MSL::Document.parse(object, &block)
    end
    
    ###
    # A MSL document.
    #
    # @!attribute [r] records
    #   @return [Array<AMDIS::MSL::Record>]
    class Document
      class << self
        ###
        # Parse MSL.  +string_or_io+ may be a String, or any object that
        # responds to _read_ and _close_ such as an IO, or StringIO.
        #
        # @param string_or_io [Object]
        # @yieldparam doc [AMDIS::MSL::Document]
        # @return [AMDIS::MSL::Document]
        def parse(string_or_io, &block)
          if string_or_io.respond_to?(:read)
            read_io(string_or_io, &block)
          else
            read_memory(string_or_io, &block)
          end
        end

        private

        ###
        # Read from object that responds to _read_ and _close_ such as an IO
        # or StringIO.
        #
        # @private
        # @param io [:read]
        # @yieldparam doc [AMDIS::MSL::Document]
        # @return [AMDIS::MSL::Document]
        def read_io(io, &block)
          read_memory(io.read, &block)
        end

        ###
        # Read from object that responds to _to_s_ such as a String.
        #
        # @private
        # @param string [:to_s]
        # @yieldparam doc [AMDIS::MSL::Document]
        # @return [AMDIS::MSL::Document]
        def read_memory(string, &block)
          records = string.to_s.scan(AMDIS::MSL::Regexen::RECORD).collect { |block_md|
            AMDIS::MSL::Record.new(*[
              block_md[0].to_i,                  # :compound_id
              block_md[1].strip,                 # :compound_name
              block_md[3].strip,                 # :molecular_formula
              block_md[4].to_f,                  # :molecular_weight
              hyphenate_cas_number(block_md[5]), # :cas_number
              block_md[6].to_f,                  # :retention_index
              block_md[7].to_f,                  # :retention_time
              block_md[8].to_f,                  # :response_factor
              block_md[9].to_f,                  # :resolution
              block_md[10].strip,                # :comment
              block_md[11].to_i,                 # :peaks_count
              block_md[12].scan(AMDIS::MSL::Regexen::PEAK).collect { |peak_md|
                AMDIS::MSL::Peak.new(*[
                  peak_md[0].to_i, # :mass_to_charge_ratio
                  peak_md[1].to_i, # :height
                ])
              },
            ])
          }
          
          doc = AMDIS::MSL::Document.new(records, &block)
          
          return doc
        end

        ###
        # Hyphenate Chemical Abstracts Service (CAS) number.
        #
        # According to the specification, a CAS number consists of three,
        # hyphen-separated segments, consisting of:
        #   * Between 2 and 6 digits;
        #   * Exactly 2 digits; and,
        #   * Exactly 1 digit.
        #
        # Since the length of the last two segments is exactly specified,
        # alignment is from right-to-left (and not left-to-right).
        #
        # This function returns nil if the CAS number is invalid (e.g., equal
        # to zero, negative, not enough digits, etc.).
        #
        # For example:
        #
        #   hyphenate_cas_number(-1) => nil
        #   hyphenate_cas_number(1118689) => '1118-68-9'
        #
        # @private
        # @param integer [:to_i]
        # @return [String, nil]
        def hyphenate_cas_number(fixnum)
          return nil unless fixnum.respond_to?(:to_i)
          
          i = fixnum.to_i
          return nil if 0 >= i
          
          if !(md = /\A([0-9]{2,6})([0-9]{2})([0-9]{1})\Z/.match(i.to_s)).nil?
            md[1..3].join('-')
          else
            nil
          end
        end
      end

      attr_reader :records

      ###
      # Default constructor.
      #
      # @param records [Array<AMDIS::MSL::Record>]
      # @yieldparam doc [AMDIS::MSL::Document]
      # @return [AMDIS::MSL::Document]
      def initialize(records = [], &block)
        @records = records
        
        if block_given?
          case block.arity
            when 1 then block.call(self)
            else        self.instance_eval(&block)
          end
        end
      end
    end

    ###
    # This module exports regular expressions for MSL documents.
    module Regexen
      ###
      # Integral number (unsigned).
      #
      # @return [Regexp]
      INT = /0|[1-9][0-9]*/.freeze

      ###
      # Floating point number (unsigned).
      #
      # @return [Regexp]
      FLOAT = /(?:#{INT})(?:#{Regexp.quote('.')}[0-9]+)?/.freeze

      ###
      # Molecular formula.
      #
      # @return [Regexp]
      MOLECULAR_FORMULA = /(?:[A-Z][a-z]*(?:#{INT})?)+/.freeze

      ###
      # MSL peak.
      #
      # @return [Regexp]
      PEAK = Regexp.new([
        Regexp.quote('('),
        '\s*',
        '(', INT, ')',
        '\s+',
        '(', INT, ')',
        '\s*',
        Regexp.quote(')'),
      ].collect(&:to_s).join).freeze

      ###
      # MSL record.
      #
      # @return [Regexp]
      RECORD = Regexp.new([
        Regexp.new([
          Regexp.quote('NAME:'),
          '\s*',
          Regexp.quote('['),
          '(', '[^', Regexp.quote(']'), ']+', ')',
          Regexp.quote(']'),
          '\s+',
          '(.+)',
          '\s+',
          Regexp.quote('['),
          '(', FLOAT, ')',
          Regexp.quote(']'),
        ].collect(&:to_s).join),
        Regexp.new([
          Regexp.quote('FORM:'),
          '\s*',
          '(', MOLECULAR_FORMULA, ')',
        ].collect(&:to_s).join),
        Regexp.new([
          Regexp.quote('MW:'),
          '\s*',
          '(', FLOAT, ')',
        ].collect(&:to_s).join),
        Regexp.new([
          Regexp.quote('CASNO:'),
          '\s*',
          '(', '(?:', INT, '|', '0+', ')', ')',
        ].collect(&:to_s).join),
        Regexp.new([
          Regexp.quote('RI:'),
          '\s*',
          '(', FLOAT, ')',
        ].collect(&:to_s).join),
        Regexp.new([
          Regexp.quote('RT:'),
          '\s*',
          '(', FLOAT, ')',
        ].collect(&:to_s).join),
        Regexp.new([
          Regexp.quote('RF:'),
          '\s*',
          '(', FLOAT, ')',
        ].collect(&:to_s).join),
        Regexp.new([
          Regexp.quote('RSN:'),
          '\s*',
          '(', FLOAT, ')',
        ].collect(&:to_s).join),
        Regexp.new([
          Regexp.quote('COMMENT:'),
          '\s*',
          '(', '.*', ')',
        ].collect(&:to_s).join),
        Regexp.new([
          Regexp.quote('Num Peaks:'),
          '\s*',
          '(', INT, ')',
          '(', '(?:\s*', PEAK, ')*', ')'
        ].collect(&:to_s).join),
      ].collect(&:to_s).join('\s+'), Regexp::MULTILINE).freeze
    end

    ###
    # A MSL record.
    #
    # @!attribute [r] compound_id
    #   @return [Fixnum]
    # @!attribute [r] compound_name
    #   @return [String]
    # @!attribute [r] molecular_formula
    #   @return [String]
    # @!attribute [r] molecular_weight
    #   @return [Float]
    # @!attribute [r] cas_number
    #   @return [String]
    # @!attribute [r] retention_index
    #   @return [Float]
    # @!attribute [r] retention_time
    #   @return [Float]
    # @!attribute [r] response_factor
    #   @return [Float]
    # @!attribute [r] resolution
    #   @return [Float]
    # @!attribute [r] comment
    #   @return [String]
    # @!attribute [r] peaks_count
    #   @return [Fixnum]
    # @!attribute [r] peaks
    #   @return [Array<AMDIS::MSL::Peak>]
    Record = Struct.new(*[
      :compound_id,
      :compound_name,
      :molecular_formula,
      :molecular_weight,
      :cas_number,
      :retention_index,
      :retention_time,
      :response_factor,
      :resolution,
      :comment,
      :peaks_count,
      :peaks,
    ])

    ###
    # A MSL peak.
    #
    # @!attribute [r] mass_to_charge_ratio
    #   @return [Fixnum]
    # @!attribute [r] height
    #   @return [Fixnum]
    Peak = Struct.new(*[
      :mass_to_charge_ratio,
      :height,
    ])
  end
end
